'use client'

import { useState } from 'react'
import { useApiKeyStore, Provider, CustomPresetKey, CustomProviderConfig } from '@/lib/stores/api-key-store-simple'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogFooter,
  DialogDescription,
} from '@/components/ui/dialog'
import {
  Key,
  Eye,
  EyeOff,
  Copy,
  Pencil,
  Trash2,
  Check,
  AlertCircle,
  Loader2
} from 'lucide-react'

// Standard provider configurations
const PROVIDER_CONFIG: Record<Exclude<Provider, 'custom'>, { name: string; description: string }> = {
  openrouter: { name: 'OpenRouter', description: 'Access multiple models' },
  anthropic: { name: 'Anthropic', description: 'Claude models' },
  openai: { name: 'OpenAI', description: 'GPT models' }
}

// Custom provider preset configurations
const CUSTOM_PRESET_CONFIG: Record<CustomPresetKey, { name: string; description: string }> = {
  deepseek: { name: 'DeepSeek', description: 'DeepSeek models' },
  qwen: { name: 'Qwen (Alibaba)', description: 'Qwen models' },
  kimi: { name: 'Kimi (Moonshot)', description: 'Moonshot models' },
  siliconflow: { name: 'SiliconFlow', description: 'SiliconFlow models' },
  minimax: { name: 'MiniMax', description: 'MiniMax models' },
  zhipuai: { name: 'Zhipu AI', description: 'GLM models' },
  manual: { name: 'Custom Endpoint', description: 'Manual configuration' }
}

function maskApiKey(key: string): string {
  if (!key || key.length < 10) return key ? '****' : ''
  return `${key.slice(0, 4)}...${key.slice(-4)}`
}

// Standard provider row component
interface StandardApiKeyRowProps {
  provider: Exclude<Provider, 'custom'>
  apiKey: string
  onEdit: () => void
  onDelete: () => void
}

function StandardApiKeyRow({ provider, apiKey, onEdit, onDelete }: StandardApiKeyRowProps) {
  const [showKey, setShowKey] = useState(false)
  const [copied, setCopied] = useState(false)
  const config = PROVIDER_CONFIG[provider]

  const handleCopy = async () => {
    if (!apiKey) return
    await navigator.clipboard.writeText(apiKey)
    setCopied(true)
    setTimeout(() => setCopied(false), 2000)
  }

  const displayValue = apiKey
    ? (showKey ? apiKey : maskApiKey(apiKey))
    : 'Not set'

  return (
    <div className="flex items-center justify-between p-3 border rounded-lg">
      <div className="flex items-center gap-3 min-w-0 flex-1">
        <Key className="h-4 w-4 text-muted-foreground flex-shrink-0" />
        <div className="min-w-0">
          <div className="font-medium text-sm">{config.name}</div>
          <div
            className={`text-xs font-mono truncate ${apiKey ? 'text-foreground' : 'text-muted-foreground'}`}
            title={showKey ? apiKey : undefined}
          >
            {displayValue}
          </div>
        </div>
      </div>

      <div className="flex items-center gap-1 flex-shrink-0">
        {apiKey && (
          <>
            <Button
              variant="ghost"
              size="icon"
              className="h-8 w-8"
              onClick={() => setShowKey(!showKey)}
              title={showKey ? 'Hide key' : 'Show key'}
            >
              {showKey ? <EyeOff className="h-4 w-4" /> : <Eye className="h-4 w-4" />}
            </Button>
            <Button
              variant="ghost"
              size="icon"
              className="h-8 w-8"
              onClick={handleCopy}
              title="Copy to clipboard"
            >
              {copied ? <Check className="h-4 w-4 text-green-500" /> : <Copy className="h-4 w-4" />}
            </Button>
          </>
        )}
        <Button
          variant="ghost"
          size="icon"
          className="h-8 w-8"
          onClick={onEdit}
          title={apiKey ? 'Edit key' : 'Add key'}
        >
          <Pencil className="h-4 w-4" />
        </Button>
        {apiKey && (
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={onDelete}
            title="Delete key"
          >
            <Trash2 className="h-4 w-4 text-red-500" />
          </Button>
        )}
      </div>
    </div>
  )
}

// Custom provider row component
interface CustomApiKeyRowProps {
  preset: CustomPresetKey
  config: CustomProviderConfig
  onEdit: () => void
  onDelete: () => void
}

function CustomApiKeyRow({ preset, config, onEdit, onDelete }: CustomApiKeyRowProps) {
  const [showKey, setShowKey] = useState(false)
  const [copied, setCopied] = useState(false)
  const presetConfig = CUSTOM_PRESET_CONFIG[preset]
  const apiKey = config.apiKey

  const handleCopy = async () => {
    if (!apiKey) return
    await navigator.clipboard.writeText(apiKey)
    setCopied(true)
    setTimeout(() => setCopied(false), 2000)
  }

  const displayValue = apiKey
    ? (showKey ? apiKey : maskApiKey(apiKey))
    : 'Not set'

  return (
    <div className="flex items-center justify-between p-3 border rounded-lg bg-muted/20">
      <div className="flex items-center gap-3 min-w-0 flex-1">
        <Key className="h-4 w-4 text-muted-foreground flex-shrink-0" />
        <div className="min-w-0">
          <div className="font-medium text-sm flex items-center gap-2">
            {presetConfig.name}
            <span className="text-xs px-1.5 py-0.5 rounded bg-primary/10 text-primary">Custom</span>
          </div>
          <div
            className={`text-xs font-mono truncate ${apiKey ? 'text-foreground' : 'text-muted-foreground'}`}
            title={showKey ? apiKey : undefined}
          >
            {displayValue}
          </div>
        </div>
      </div>

      <div className="flex items-center gap-1 flex-shrink-0">
        {apiKey && (
          <>
            <Button
              variant="ghost"
              size="icon"
              className="h-8 w-8"
              onClick={() => setShowKey(!showKey)}
              title={showKey ? 'Hide key' : 'Show key'}
            >
              {showKey ? <EyeOff className="h-4 w-4" /> : <Eye className="h-4 w-4" />}
            </Button>
            <Button
              variant="ghost"
              size="icon"
              className="h-8 w-8"
              onClick={handleCopy}
              title="Copy to clipboard"
            >
              {copied ? <Check className="h-4 w-4 text-green-500" /> : <Copy className="h-4 w-4" />}
            </Button>
          </>
        )}
        <Button
          variant="ghost"
          size="icon"
          className="h-8 w-8"
          onClick={onEdit}
          title={apiKey ? 'Edit key' : 'Add key'}
        >
          <Pencil className="h-4 w-4" />
        </Button>
        {apiKey && (
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={onDelete}
            title="Delete key"
          >
            <Trash2 className="h-4 w-4 text-red-500" />
          </Button>
        )}
      </div>
    </div>
  )
}

type EditTarget =
  | { type: 'standard'; provider: Exclude<Provider, 'custom'> }
  | { type: 'custom'; preset: CustomPresetKey }
  | null

type DeleteTarget = EditTarget

export function ApiKeyManager() {
  const {
    apiKeys,
    setApiKey,
    clearApiKey,
    customProviders,
    setCustomProviderKey,
    clearCustomProviderKey
  } = useApiKeyStore()

  const [editTarget, setEditTarget] = useState<EditTarget>(null)
  const [editValue, setEditValue] = useState('')
  const [deleteTarget, setDeleteTarget] = useState<DeleteTarget>(null)
  const [isSaving, setIsSaving] = useState(false)

  const standardProviders: Array<Exclude<Provider, 'custom'>> = ['openrouter', 'anthropic', 'openai']

  // Get custom providers that have keys configured
  const configuredCustomProviders = (Object.entries(customProviders) as Array<[CustomPresetKey, CustomProviderConfig]>)
    .filter(([_, config]) => config.apiKey)

  const handleEditStandard = (provider: Exclude<Provider, 'custom'>) => {
    setEditTarget({ type: 'standard', provider })
    setEditValue(apiKeys[provider] || '')
  }

  const handleEditCustom = (preset: CustomPresetKey) => {
    setEditTarget({ type: 'custom', preset })
    setEditValue(customProviders[preset]?.apiKey || '')
  }

  const handleSave = async () => {
    if (!editTarget) return
    setIsSaving(true)
    try {
      if (editTarget.type === 'standard') {
        await setApiKey(editValue, editTarget.provider)
      } else {
        await setCustomProviderKey(editTarget.preset, editValue)
      }
      setEditTarget(null)
      setEditValue('')
    } finally {
      setIsSaving(false)
    }
  }

  const handleDelete = async () => {
    if (!deleteTarget) return
    setIsSaving(true)
    try {
      if (deleteTarget.type === 'standard') {
        await clearApiKey(deleteTarget.provider)
      } else {
        await clearCustomProviderKey(deleteTarget.preset)
      }
      setDeleteTarget(null)
    } finally {
      setIsSaving(false)
    }
  }

  const handleEditClose = (open: boolean) => {
    if (!open) {
      setEditTarget(null)
      setEditValue('')
    }
  }

  const handleDeleteClose = (open: boolean) => {
    if (!open) {
      setDeleteTarget(null)
    }
  }

  const getEditTitle = () => {
    if (!editTarget) return ''
    if (editTarget.type === 'standard') {
      const hasKey = apiKeys[editTarget.provider]
      return `${hasKey ? 'Edit' : 'Add'} ${PROVIDER_CONFIG[editTarget.provider].name} API Key`
    } else {
      const hasKey = customProviders[editTarget.preset]?.apiKey
      return `${hasKey ? 'Edit' : 'Add'} ${CUSTOM_PRESET_CONFIG[editTarget.preset].name} API Key`
    }
  }

  const getDeleteName = () => {
    if (!deleteTarget) return ''
    if (deleteTarget.type === 'standard') {
      return PROVIDER_CONFIG[deleteTarget.provider].name
    } else {
      return CUSTOM_PRESET_CONFIG[deleteTarget.preset].name
    }
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Key className="h-5 w-5" />
          API Keys
        </CardTitle>
      </CardHeader>
      <CardContent className="space-y-3">
        {/* Standard Providers */}
        {standardProviders.map((provider) => (
          <StandardApiKeyRow
            key={provider}
            provider={provider}
            apiKey={apiKeys[provider]}
            onEdit={() => handleEditStandard(provider)}
            onDelete={() => setDeleteTarget({ type: 'standard', provider })}
          />
        ))}

        {/* Custom Providers with configured keys */}
        {configuredCustomProviders.length > 0 && (
          <>
            <div className="text-xs text-muted-foreground pt-2 border-t mt-4">
              Custom Providers
            </div>
            {configuredCustomProviders.map(([preset, config]) => (
              <CustomApiKeyRow
                key={preset}
                preset={preset}
                config={config}
                onEdit={() => handleEditCustom(preset)}
                onDelete={() => setDeleteTarget({ type: 'custom', preset })}
              />
            ))}
          </>
        )}

        {/* Edit Dialog */}
        <Dialog open={!!editTarget} onOpenChange={handleEditClose}>
          <DialogContent>
            <DialogHeader>
              <DialogTitle>{getEditTitle()}</DialogTitle>
              <DialogDescription>
                Enter your API key. It will be stored securely.
                {editTarget?.type === 'custom' && ' (Custom provider keys are synced to your account.)'}
              </DialogDescription>
            </DialogHeader>
            <div className="py-4">
              <Input
                type="password"
                placeholder="Enter API key..."
                value={editValue}
                onChange={(e) => setEditValue(e.target.value)}
                autoFocus
              />
            </div>
            <DialogFooter>
              <Button variant="outline" onClick={() => setEditTarget(null)}>
                Cancel
              </Button>
              <Button onClick={handleSave} disabled={isSaving || !editValue.trim()}>
                {isSaving && <Loader2 className="h-4 w-4 animate-spin mr-2" />}
                Save
              </Button>
            </DialogFooter>
          </DialogContent>
        </Dialog>

        {/* Delete Confirmation Dialog */}
        <Dialog open={!!deleteTarget} onOpenChange={handleDeleteClose}>
          <DialogContent>
            <DialogHeader>
              <DialogTitle className="flex items-center gap-2">
                <AlertCircle className="h-5 w-5 text-red-500" />
                Delete API Key
              </DialogTitle>
              <DialogDescription>
                Are you sure you want to delete the {getDeleteName()} API key?
                This action cannot be undone.
              </DialogDescription>
            </DialogHeader>
            <DialogFooter>
              <Button variant="outline" onClick={() => setDeleteTarget(null)}>
                Cancel
              </Button>
              <Button variant="destructive" onClick={handleDelete} disabled={isSaving}>
                {isSaving && <Loader2 className="h-4 w-4 animate-spin mr-2" />}
                Delete
              </Button>
            </DialogFooter>
          </DialogContent>
        </Dialog>
      </CardContent>
    </Card>
  )
}
