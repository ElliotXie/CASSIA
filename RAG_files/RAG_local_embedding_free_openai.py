
#%%capture --no-stderr
#%pip install -U langchain_community tiktoken langchainhub scikit-learn langchain langgraph tavily-python  nomic[local] langchain-nomic langchain_openai
#
# %%

import getpass
import os
import pandas as pd
from typing import List
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.document_loaders import DirectoryLoader, CSVLoader, UnstructuredExcelLoader
from langchain_community.vectorstores import SKLearnVectorStore
from langchain.schema import Document
from langchain.prompts import PromptTemplate
from langchain_community.chat_models import ChatOllama
from langchain_core.output_parsers import StrOutputParser
from langchain_community.embeddings import HuggingFaceEmbeddings

def set_env(key: str):
    os.environ[key] = getpass.getpass(f"{key}:")

local_llm = "llama3.1:latest"

# Function to load CSV and Excel files
def load_csv_excel_files(folder_path: str) -> List[Document]:
    documents = []
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if filename.endswith('.csv'):
            loader = CSVLoader(file_path)
            documents.extend(loader.load())
        elif filename.endswith(('.xlsx', '.xls')):
            loader = UnstructuredExcelLoader(file_path)
            documents.extend(loader.load())
    return documents

# Specify the path to your local folder containing papers and data files
local_folder_path = "/Users/lcheng74/RAG_Papers"


text_docs = DirectoryLoader(local_folder_path, glob="**/*.txt").load()
pdf_docs = DirectoryLoader(local_folder_path, glob="**/*.pdf").load()

data_docs = load_csv_excel_files(local_folder_path)

# Combine all documents
all_docs = text_docs + data_docs + pdf_docs

docs_list = [doc for doc in all_docs if isinstance(doc, Document)]

text_splitter = RecursiveCharacterTextSplitter.from_tiktoken_encoder(
    chunk_size=250, chunk_overlap=0
)

doc_splits = text_splitter.split_documents(docs_list)

embedding = HuggingFaceEmbeddings(model_name="all-MiniLM-L6-v2")

vectorstore = SKLearnVectorStore.from_documents(
    documents=doc_splits,
    embedding=embedding,
)

retriever = vectorstore.as_retriever(k=4)

llm = ChatOllama(model=local_llm, temperature=0)

prompt = PromptTemplate(
    template="""You are an expert in cell biology and tissue analysis. Your task is to help find markers and key features of different cell types in various tissues, optionally under specific conditions.

Cell Type: {cell_type}
Tissue Type: {tissue_type}
Condition: {condition}

Using the following documents, identify and summarize the key markers and features of the specified cell type in the given tissue. If a condition is specified, focus on how it affects the cell type's characteristics, key works might include: overprepresent, underrepresent, differential expressed, highly expressed, low, upregulated, downregulated, etc.

Documents: {documents}

Provide a concise summary of:
1. Specific markers (e.g., genes, proteins) associated with this cell type in the given tissue
2. Key features or functions of this cell type in the tissue
3. If a condition is specified, how it impacts the cell type's markers or features

Answer:
""",
    input_variables=["cell_type", "tissue_type", "condition", "documents"],
)

rag_chain = prompt | llm | StrOutputParser()


def get_rag_response(cell_type: str, tissue_type: str, question: str, condition: str = ""):
    # Combine cell type and tissue type for retrieval
    retrieval_query = f"{cell_type} {tissue_type}"
    if condition:
        retrieval_query += f" {condition}"
    
    docs = retriever.invoke(retrieval_query)
    
    return rag_chain.invoke({
        "cell_type": cell_type,
        "tissue_type": tissue_type,
        "condition": condition if condition else "No specific condition specified",
        "documents": docs
    })
    


cell_type = "SSC"
tissue_type = "colon"
question = "What are the markers and key features of SSC in the colon?"
condition = ""  # Optional, can be left empty

response = get_rag_response(cell_type, tissue_type, question, condition)
print(response)


# %%
