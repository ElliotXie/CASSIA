<div align="center">

<img src="CASSIA_python/logo2.png" width="200" style="vertical-align: middle;" />

[English](README.md) | [中文](README_CN.md)

</div>

**CASSIA** 是一个基于multi-agent多智能体的大型语言模型工具，用于快速，准确，简单地进行单细胞的可解释分群注释。

🌐 [体验 CASSIA 网页界面](https://cassiacell.com/) - 仅提供最基础的CASSIA功能

📚 [完整 R 文档及示例/英文版/最新](https://cassia-documentation-en-new.vercel.app/)

📚 [完整 R 文档/中文版](https://cassia-documentation-cn.vercel.app/)

📝 [Python 工作流示例](https://github.com/ElliotXie/CASSIA/blob/main/CASSIA_example/CASSIA_python_tutorial.ipynb)

🤖 [模型注释能力排行榜](https://sc-llm-benchmark.com/methods/cassia)



## 📰 更新

> **2025-05-05**  
> 📊 **CASSIA注释基准测试平台现已上线！**  
> 本次更新推出了一个全新的基准测试平台，用于评估不同大型语言模型在单细胞注释任务中的表现与成本。  
> **LLaMA4 Maverick、Gemini 2.5 Flash 和 DeepSeekV3 是目前性能与成本最均衡的模型，且几乎免费使用！**  
> 🔧 新增"自动合并功能"，可统一输出不同层级的细胞类型标签，大幅简化子聚类分析流程。  
> 🐛 修复了注释增强代理中的一个错误，提高了低质量注释的优化效果。

> **2025-04-19**  
> 🔄 **CASSIA 添加了重试机制和优化的报告存储！**  
> 最新更新引入了失败任务的自动重试机制，并优化了报告的存储方式，使访问和管理更加便捷。  
> 🎨 **完成CASSIA标志设计！**

> **2025-04-17**  
> 🚀 **CASSIA 现在支持自动单细胞注释基准测试！**  
> 最新更新引入了一个新功能，可以实现完全自动化的单细胞注释基准测试。结果由 LLMs 自动评估，性能与人类专家相当。  
> **专门的基准测试网站即将推出—敬请期待！**


## 🏗️ 安装 (R 语言，Python教程请访问[这里](https://github.com/ElliotXie/CASSIA/blob/main/CASSIA_example/CASSIA_python_tutorial.ipynb))

```R
# 安装依赖
install.packages("devtools")
install.packages("reticulate")

# 安装 CASSIA
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")
```

***注意：如果环境第一次没有正确设置，请重启R并运行以下代码***

```R
library(CASSIA)
setup_cassia_env()
```

### 🔑 设置 API 密钥

获取API密钥大约需要3分钟时间。

针对国内用户，我们强烈推荐使用 OpenRouter ，因为OpenAI和Anthropic都对国内访问有限制，使用OpenRouter可以通过单个 API 密钥访问大多数模型。

```R
# 对于 OpenRouter
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)

# 对于 OpenAI
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)

# 对于 Anthropic
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)
```


- **API 提供商指南:**

	- [如何获取 OpenRouter API 密钥](https://openrouter.ai/settings/keys)
 	- [如何充值OpenRouter](https://zhuanlan.zhihu.com/p/1898753591528908109)



## 🧬 示例数据

CASSIA 包含两种格式的示例标记数据：
```R
# 加载示例数据
markers_unprocessed <- loadExampleMarkers(processed = FALSE)  # Seurat findallmarkers 输出文件
markers_processed <- loadExampleMarkers(processed = TRUE)     # 处理后格式
```

## ⚙️ 流程使用

```R
# 默认提供商设置为OpenRouter

runCASSIA_pipeline(
    output_file_name = "cassia_test",            # Base name for output files
    tissue = "Large Intestine",                   # Tissue type (e.g., "brain")
    species = "Human",              		 # Species (e.g., "human")
    marker = "markers_unprocessed",               # Marker data from findallmarker
    max_workers = 4                              # Number of parallel workers
)
```

## 🤖 支持的模型

您可以为注释和评分选择任何模型。下面列出了一些经典模型。OpenRouter支持当前大多数流行的模型，尽管有些模型在CASSIA论文中尚未进行全面基准测试 — 欢迎尝试。


### OpenRouter
- `google/gemini-2.5-flash-preview`: 最好的低费率大模型之一（最推荐）
- `deepseek/deepseek-chat-v3-0324`: 最好的开源大模型之一，经常给出非常详细的注释（推荐）
- `deepseek/deepseek-chat-v3-0324:free`: 免费但速度较慢

### OpenAI
- `gpt-4o`: 用于文章的基准测试

### Anthropic
- `claude-3-7-sonnet-latest`: 最新的高性能模型

## 📤 输出

流程生成四个关键文件：
1. 完整注释结果表格
2. 注释摘要网页报告
3. 注释增强智能体报告

## 🧰 故障排除

### 身份验证（错误 401）
```R
# 检查 API 密钥是否正确设置
key <- Sys.getenv("ANTHROPIC_API_KEY")
print(key)  # 输出结果不应该为空

# 如果需要，重置 API 密钥
setLLMApiKey("your_api_key", provider = "anthropic", persist = TRUE)
```

### 文件错误
- 必要时使用绝对路径
- 检查文件权限
- 确保文件未在其他程序中打开

### 最佳实践
- 保持API密钥安全
- 维持足够的API额度


注意：此 README 仅涵盖了基本的 CASSIA 功能。有关包括高级功能和详细示例在内的完整教程，请访问：
[CASSIA 完整教程](https://cassia-documentation-en-new.vercel.app/)。

## 📖 引用

📖 [阅读我们的预印本 (v2, 最新版)](https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2)
 
📖 [原始预印本 (v1, 历史版本)](https://www.biorxiv.org/content/10.1101/2024.12.04.626476v1)

CASSIA: a multi-agent large language model for reference-free, interpretable, and automated cell annotation of single-cell RNA-sequencing data  
Elliot Xie, Lingxin Cheng, Jack Shireman, Yujia Cai, Jihua Liu, Chitrasen Mohanty, Mahua Dey, Christina Kendziorski  
bioRxiv 2024.12.04.626476; doi: https://doi.org/10.1101/2024.12.04.626476

## 📬 联系方式

如有任何问题或需要帮助，欢迎随时邮件联系我们，我们一定会尽力协助：
**xie227@wisc.edu** 
如果您觉得我们的项目对您有帮助，请留下一个⭐或者将项目分享给你的朋友，感激不尽！
