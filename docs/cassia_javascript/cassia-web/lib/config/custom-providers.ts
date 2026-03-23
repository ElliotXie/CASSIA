/**
 * Custom provider presets for third-party OpenAI-compatible API endpoints.
 * Shared across ApiKeyInput, AgentModelSelector, and page.tsx.
 *
 * These contain API base URLs and provider-specific model lists,
 * NOT managed by model_settings.json (which handles standard providers only).
 */

export const CUSTOM_PROVIDER_PRESETS = {
  deepseek: {
    name: 'DeepSeek',
    baseUrl: 'https://api.deepseek.com',
    models: ['deepseek-chat', 'deepseek-reasoner'],
    helpUrl: 'https://platform.deepseek.com/api_keys'
  },
  qwen: {
    name: 'Qwen (Alibaba)',
    baseUrl: 'https://dashscope-intl.aliyuncs.com/compatible-mode/v1',
    models: ['qwen-max', 'qwen-plus', 'qwen-turbo'],
    helpUrl: 'https://dashscope.console.aliyun.com/apiKey'
  },
  kimi: {
    name: 'Kimi (Moonshot)',
    baseUrl: 'https://api.moonshot.cn/v1',
    models: ['kimi-k2.5'],
    helpUrl: 'https://platform.moonshot.cn/console/api-keys'
  },
  siliconflow: {
    name: 'SiliconFlow',
    baseUrl: 'https://api.siliconflow.cn/v1',
    models: ['deepseek-ai/DeepSeek-V3.2'],
    helpUrl: 'https://cloud.siliconflow.cn/account/ak'
  },
  minimax: {
    name: 'MiniMax',
    baseUrl: 'https://api.minimax.io/v1',
    models: ['MiniMax-M2'],
    helpUrl: 'https://platform.minimaxi.com/user-center/basic-information/interface-key'
  },
  zhipuai: {
    name: 'Zhipu AI',
    baseUrl: 'https://open.bigmodel.cn/api/paas/v4',
    models: ['glm-5'],
    helpUrl: 'https://open.bigmodel.cn/usercenter/apikeys'
  },
  manual: {
    name: 'Manual Entry',
    baseUrl: '',
    models: [],
    helpUrl: '#'
  }
} as const;

export type CustomPresetKey = keyof typeof CUSTOM_PROVIDER_PRESETS;
