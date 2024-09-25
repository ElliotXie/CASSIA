# %%
import getpass
import os
from typing import List
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.document_loaders import DirectoryLoader, CSVLoader, UnstructuredExcelLoader
from langchain.schema import Document
from langchain.prompts import PromptTemplate
from langchain_community.chat_models import ChatOllama
from langchain_core.output_parsers import StrOutputParser
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain.retrievers import BM25Retriever
from langchain_community.vectorstores import FAISS
import networkx as nx
from sentence_transformers import util
import torch

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

# Load text documents
text_docs = DirectoryLoader(local_folder_path, glob="**/*.pdf").load()

# Load CSV and Excel files
data_docs = load_csv_excel_files(local_folder_path)

# Combine all documents
all_docs = text_docs + data_docs

# Ensure all docs are Document objects
docs_list = [doc for doc in all_docs if isinstance(doc, Document)]

# Initialize a text splitter with specified chunk size and overlap
text_splitter = RecursiveCharacterTextSplitter.from_tiktoken_encoder(
    chunk_size=250, chunk_overlap=50
)

# Split the documents into chunks
doc_splits = text_splitter.split_documents(docs_list)

# Use SentenceTransformers embedding
embedding_model = HuggingFaceEmbeddings(model_name="all-MiniLM-L6-v2")

# Create FAISS vector store
vectorstore = FAISS.from_documents(doc_splits, embedding_model)

# Create BM25 retriever for keyword search
bm25_retriever = BM25Retriever.from_documents(doc_splits)
bm25_retriever.k = 4  # Retrieve top 4 documents

# Function to create graph from document chunks
def create_graph(chunks: List[Document]) -> nx.Graph:
    G = nx.Graph()
    embeddings = embedding_model.embed_documents([chunk.page_content for chunk in chunks])
    
    for i, chunk in enumerate(chunks):
        G.add_node(i, text=chunk.page_content, embedding=embeddings[i])
    
    # Create edges based on cosine similarity
    for i in range(len(chunks)):
        for j in range(i+1, len(chunks)):
            similarity = util.cos_sim(embeddings[i], embeddings[j]).item()
            if similarity > 0.5:  # Adjust threshold as needed
                G.add_edge(i, j, weight=similarity)
    
    return G

# Create graph
graph = create_graph(doc_splits)

# Hybrid retrieval function
def hybrid_retrieve(query: str, k: int = 4) -> List[str]:
    # Vector search
    vector_results = vectorstore.similarity_search(query, k=k)
    
    # BM25 search (updated to use invoke)
    bm25_results = bm25_retriever.invoke(query)
    
    # Graph-based search
    query_embedding = embedding_model.embed_query(query)
    graph_similarities = []
    for node in graph.nodes():
        node_embedding = graph.nodes[node]['embedding']
        similarity = util.cos_sim(query_embedding, node_embedding).item()
        graph_similarities.append((node, similarity))
    
    top_graph_nodes = sorted(graph_similarities, key=lambda x: x[1], reverse=True)[:k]
    graph_results = [Document(page_content=graph.nodes[node]['text']) for node, _ in top_graph_nodes]
    
    # Combine results
    all_results = vector_results + bm25_results + graph_results
    
    # Remove duplicates while preserving order
    seen = set()
    unique_results = []
    for doc in all_results:
        if doc.page_content not in seen:
            seen.add(doc.page_content)
            unique_results.append(doc.page_content)
    
    return unique_results[:k]

llm = ChatOllama(model=local_llm, temperature=0)

prompt = PromptTemplate(
    template="""You are an expert in cell biology and tissue analysis. Your task is to help find markers and key features of different cell types in various tissues, optionally under specific conditions.

Cell Type: {cell_type}
Tissue Type: {tissue_type}
Condition: {condition}

Using the following documents, identify and summarize the key markers and features of the specified cell type in the given tissue. If a condition is specified, focus on how it affects the cell type's characteristics, key words might include: overpresent, underpresent, differentially expressed, highly expressed, low, upregulated, downregulated, etc.

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

def get_hybrid_rag_response(cell_type: str, tissue_type: str, question: str, condition: str = ""):
    # Combine cell type and tissue type for retrieval
    retrieval_query = f"{cell_type} {tissue_type}"
    if condition:
        retrieval_query += f" {condition}"
    
    # Use hybrid retrieval to get relevant documents
    docs = hybrid_retrieve(retrieval_query)
    
    return rag_chain.invoke({
        "cell_type": cell_type,
        "tissue_type": tissue_type,
        "condition": condition if condition else "No specific condition specified",
        "documents": docs
    })

# Example usage
if __name__ == "__main__":
    cell_type = "SSC"
    tissue_type = "colon"
    question = "What are the markers and key features of SSC in the colon?"
    condition = ""  # Optional, can be left empty

    response = get_hybrid_rag_response(cell_type, tissue_type, question, condition)
    print(response)
# %%
