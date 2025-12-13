'use client'

import { useState } from 'react'
import { useApiKeyStore, Provider } from '@/lib/stores/api-key-store-simple'
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

const PROVIDER_CONFIG: Record<Provider, { name: string; description: string }> = {
  openrouter: { name: 'OpenRouter', description: 'Access multiple models' },
  anthropic: { name: 'Anthropic', description: 'Claude models' },
  openai: { name: 'OpenAI', description: 'GPT models' },
  custom: { name: 'Custom', description: 'Local only' }
}

function maskApiKey(key: string): string {
  if (!key || key.length < 10) return key ? '****' : ''
  return `${key.slice(0, 4)}...${key.slice(-4)}`
}

interface ApiKeyRowProps {
  provider: Provider
  apiKey: string
  onEdit: (provider: Provider) => void
  onDelete: (provider: Provider) => void
}

function ApiKeyRow({ provider, apiKey, onEdit, onDelete }: ApiKeyRowProps) {
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
          <div className="font-medium text-sm flex items-center gap-2">
            {config.name}
            {provider === 'custom' && (
              <span className="text-xs text-muted-foreground">(local only)</span>
            )}
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
          onClick={() => onEdit(provider)}
          title={apiKey ? 'Edit key' : 'Add key'}
        >
          <Pencil className="h-4 w-4" />
        </Button>
        {apiKey && (
          <Button
            variant="ghost"
            size="icon"
            className="h-8 w-8"
            onClick={() => onDelete(provider)}
            title="Delete key"
          >
            <Trash2 className="h-4 w-4 text-red-500" />
          </Button>
        )}
      </div>
    </div>
  )
}

export function ApiKeyManager() {
  const { apiKeys, setApiKey, clearApiKey, isLoading } = useApiKeyStore()

  const [editingProvider, setEditingProvider] = useState<Provider | null>(null)
  const [editValue, setEditValue] = useState('')
  const [deleteConfirm, setDeleteConfirm] = useState<Provider | null>(null)
  const [isSaving, setIsSaving] = useState(false)

  const providers: Provider[] = ['openrouter', 'anthropic', 'openai', 'custom']

  const handleEdit = (provider: Provider) => {
    setEditingProvider(provider)
    setEditValue(apiKeys[provider] || '')
  }

  const handleSave = async () => {
    if (!editingProvider) return
    setIsSaving(true)
    try {
      await setApiKey(editValue, editingProvider)
      setEditingProvider(null)
      setEditValue('')
    } finally {
      setIsSaving(false)
    }
  }

  const handleDelete = async () => {
    if (!deleteConfirm) return
    setIsSaving(true)
    try {
      await clearApiKey(deleteConfirm)
      setDeleteConfirm(null)
    } finally {
      setIsSaving(false)
    }
  }

  const handleEditClose = (open: boolean) => {
    if (!open) {
      setEditingProvider(null)
      setEditValue('')
    }
  }

  const handleDeleteClose = (open: boolean) => {
    if (!open) {
      setDeleteConfirm(null)
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
        {providers.map((provider) => (
          <ApiKeyRow
            key={provider}
            provider={provider}
            apiKey={apiKeys[provider]}
            onEdit={handleEdit}
            onDelete={(p) => setDeleteConfirm(p)}
          />
        ))}

        {/* Edit Dialog */}
        <Dialog open={!!editingProvider} onOpenChange={handleEditClose}>
          <DialogContent>
            <DialogHeader>
              <DialogTitle>
                {editingProvider && apiKeys[editingProvider] ? 'Edit' : 'Add'}{' '}
                {editingProvider && PROVIDER_CONFIG[editingProvider]?.name} API Key
              </DialogTitle>
              <DialogDescription>
                Enter your API key. It will be stored securely.
                {editingProvider === 'custom' && ' (Custom keys are stored locally only.)'}
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
              <Button variant="outline" onClick={() => setEditingProvider(null)}>
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
        <Dialog open={!!deleteConfirm} onOpenChange={handleDeleteClose}>
          <DialogContent>
            <DialogHeader>
              <DialogTitle className="flex items-center gap-2">
                <AlertCircle className="h-5 w-5 text-red-500" />
                Delete API Key
              </DialogTitle>
              <DialogDescription>
                Are you sure you want to delete the{' '}
                {deleteConfirm && PROVIDER_CONFIG[deleteConfirm]?.name} API key?
                This action cannot be undone.
              </DialogDescription>
            </DialogHeader>
            <DialogFooter>
              <Button variant="outline" onClick={() => setDeleteConfirm(null)}>
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
