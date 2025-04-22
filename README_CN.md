# CASSIA <img src="CASSIA_python/logo2.png" align="right" width="150" style="vertical-align: middle;" />

[English](README.md) | [中文](README_CN.md)

CASSIA (用于单细胞可解释注释的协作智能体系统) 是一种利用多智能体大型语言模型（LLMs）增强细胞类型注释的工具。

🌐 [体验 CASSIA 网页界面](https://cassiacell.com/) - 基础 CASSIA 功能的网络界面

📝 [R 工作流示例](https://github.com/ElliotXie/CASSIA/blob/main/CASSIA_example/CASSIA_tutorial_final.Rmd)

📚 [完整 R 文档](https://cassia-true-final-4.vercel.app/)

📝 [Python 工作流示例](https://github.com/ElliotXie/CASSIA/blob/main/CASSIA_example/CASSIA_python_tutorial.ipynb)

📖 [最新预印本 (v2, 最新版)](https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2)
 
📖 [原始预印本 (v1, 历史版本)](https://www.biorxiv.org/content/10.1101/2024.12.04.626476v1)


## 📰 更新

> **2025-04-19**  
> 🔄 **CASSIA 添加了重试机制和优化的报告存储！**  
> 最新更新引入了失败任务的自动重试机制，并优化了报告的存储方式，使访问和管理更加便捷。  
> 🎨 **CASSIA 标志已经设计并添加到项目中！**

> **2025-04-17**  
> 🚀 **CASSIA 现在支持自动单细胞注释基准测试！**  
> 最新更新引入了一个新功能，可以实现完全自动化的单细胞注释基准测试。结果由 LLMs 自动评估，性能与人类专家相当。  
> **专门的基准测试网站即将推出—敬请期待！**


## 🏗️ 安装 (R)

从 GitHub 安装
```R
# 安装依赖
install.packages("devtools")
install.packages("reticulate")

# 安装 CASSIA
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")
```

### 🔑 设置 API 密钥

我们建议从 OpenRouter 开始，因为它可以通过单个 API 密钥访问大多数模型。虽然价格略贵且偶尔不稳定，但它提供了更大的便利性。对于生产用途，通过 OpenAI 或 Anthropic 直接访问提供了更好的稳定性。

请注意，在某些国家，OpenAI 和 Anthropic 可能被禁止。在这些情况下，用户可以使用 OpenRouter 代替。

```R
# 对于 OpenAI
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)

# 对于 Anthropic
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)

# 对于 OpenRouter
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)
```


- **API 提供商指南:**
	- [如何获取 OpenAI API 密钥](https://platform.openai.com/api-keys)
	- [如何获取 Anthropic API 密钥](https://console.anthropic.com/settings/keys)
	- [如何获取 OpenRouter API 密钥](https://openrouter.ai/settings/keys)
    - [OpenAI API 文档](https://beta.openai.com/docs/)
    - [Anthropic API 文档](https://docs.anthropic.com/)
    - [OpenRouter API 文档](https://openrouter.ai/docs/quick-start)


## 🧬 示例数据

CASSIA 包含两种格式的示例标记数据：
```R
# 加载示例数据
markers_unprocessed <- loadExampleMarkers(processed = FALSE)  # 直接 Seurat 输出
markers_processed <- loadExampleMarkers(processed = TRUE)     # 处理后格式
```

## ⚙️ 流程使用

```R
runCASSIA_pipeline(
    output_file_name,     # 输出文件的基本名称
    tissue,               # 组织类型（例如，"brain"）
    species,              # 物种（例如，"human"）
    marker,               # 来自 findallmarker 的标记数据
    max_workers = 4,      # 并行工作者数量
    annotation_model = "gpt-4o",                    # 注释模型
    annotation_provider = "openai",                 # 注释提供商
    score_model = "anthropic/claude-3.5-sonnet",    # 评分模型
    score_provider = "openrouter",                  # 评分提供商
    annotationboost_model="anthropic/claude-3.5-sonnet", # 注释增强模型
    annotationboost_provider="openrouter", # 注释增强提供商
    score_threshold = 75,                          # 最低可接受分数
    additional_info = NULL                         # 可选上下文信息
)
```

## 🤖 支持的模型

您可以为注释和评分选择任何模型。下面列出了一些经典模型。大多数当前流行的模型都被 OpenRouter 支持，尽管它们还没有在 CASSIA 论文中进行广泛的基准测试——随时尝试它们。

### OpenAI（最常见）
- `gpt-4o`（推荐）：性能和成本平衡
- `o1-mini`：高级推理能力（成本更高）

### Anthropic
- `claude-3-5-sonnet-20241022`：高性能模型
- `claude-3-7-sonnet-latest`：最新模型

### OpenRouter
- `anthropic/claude-3.5-sonnet`：高访问限制的 Claude 访问
- `openai/gpt-4o-2024-11-20`：GPT-4o 的替代访问
- `meta-llama/llama-3.2-90b-vision-instruct`：经济实惠的开源选项
- `deepseek/deepseek-chat-v3-0324`：非常经济实惠且与 GPT-4o 相当

## 📤 输出

流程生成四个关键文件：
1. 初始注释结果
2. 带推理的质量评分
3. 摘要报告
4. 注释增强报告

## 🧰 故障排除

### 身份验证（错误 401）
```R
# 检查 API 密钥是否正确设置
key <- Sys.getenv("ANTHROPIC_API_KEY")
print(key)  # 不应该为空

# 如果需要，重置 API 密钥
setLLMApiKey("your_api_key", provider = "anthropic", persist = TRUE)
```

### 文件错误
- 必要时使用绝对路径
- 检查文件权限
- 确保文件未在其他程序中打开
- 验证磁盘空间是否充足

### 最佳实践
- 保持 API 密钥安全
- 维持足够的 API 积分
- 在覆盖文件之前备份数据
- 仔细检查文件路径和权限

注意：此 README 涵盖了基本的 CASSIA 功能。有关包括高级功能和详细示例在内的完整教程，请访问：
[CASSIA 完整教程](https://cassia-true-final-4.vercel.app/)。

## 📖 引用

CASSIA: a multi-agent large language model for reference free, interpretable, and automated cell annotation of single-cell RNA-sequencing data  
Elliot Xie, Lingxin Cheng, Jack Shireman, Yujia Cai, Jihua Liu, Chitrasen Mohanty, Mahua Dey, Christina Kendziorski  
bioRxiv 2024.12.04.626476; doi: https://doi.org/10.1101/2024.12.04.626476 