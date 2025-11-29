---
title: 故障排除
---

这里我们列出了用户可能遇到的一些常见错误。如果您有任何其他问题，请随时在最后一部分发表评论！

## 认证错误 (Error 401)

### 常见原因
* API 密钥未正确加载
* API 密钥无效或过期
* API 额度不足

### 解决方案
```python
import os
# 验证 API 密钥是否正确设置
key = os.environ.get("ANTHROPIC_API_KEY")
print(key)  # 应显示您的 API 密钥，而不是 None 或空字符串

# 如果需要，重置 API 密钥
import CASSIA
CASSIA.set_api_key("your_anthropic_api_key", provider = "anthropic")
```

确保：
* 验证您的平台账户上有足够的额度
* 如有必要，通过平台的计费门户充值

## 文件未找到错误

### 常见原因
* 文件路径不正确
* 缺少输入文件
* CASSIA 自动文件匹配问题

### 解决方案
```python
import os
# 验证文件路径正确且文件存在
os.path.exists("your_input_file.csv")  # 如果文件存在则返回 True
```

指定路径时：
* 必要时使用绝对路径：
  * 不正确: "data/input.csv"
  * 正确: "C:/Users/YourName/Project/data/input.csv"
* 确保输入和输出文件名符合 CASSIA 的预期：
  * 严格遵循命名约定
  * 检查文件扩展名是否符合所需格式

## 权限错误

### 常见原因
* 文件已存在且被锁定/正在使用
* 写入权限不足
* 目录访问限制

### 解决方案
```python
import os
# 如果文件存在，在继续之前将其删除
if os.path.exists("output.csv"):
    os.remove("output.csv")
```

其他步骤：
* 检查文件是否已在另一个程序中打开
* 关闭任何使用该文件的程序
* 等待几秒钟后重试
* 如果需要，使用不同的输出文件名

## 最佳实践
* 妥善保管 API 密钥，切勿分享
* 在覆盖文件之前务必备份重要数据
* 在运行操作之前仔细检查文件路径和权限
* 保持足够的 API 额度以进行不间断的处理

