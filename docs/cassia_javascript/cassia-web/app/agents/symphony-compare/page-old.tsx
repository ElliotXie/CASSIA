"use client";

import { useState, ChangeEvent, FC, ReactNode, useEffect } from 'react';
import { symphonyCompare } from '@/lib/cassia/symphonyCompare';
import { useApiKeyStore } from '@/lib/stores/api-key-store';

// Basic UI components with types
interface ButtonProps {
    children: ReactNode;
    onClick: () => void;
    disabled: boolean;
}

const Button: FC<ButtonProps> = ({ children, onClick, disabled }) => (
  <button onClick={onClick} disabled={disabled} style={{ padding: '10px 20px', margin: '5px', cursor: disabled ? 'not-allowed' : 'pointer' }}>
    {children}
  </button>
);

interface InputProps {
    type: string;
    value: string | number;
    onChange: (e: ChangeEvent<HTMLInputElement>) => void;
    placeholder?: string;
    className?: string;
    min?: string;
    max?: string;
    step?: string;
}

const Input: FC<InputProps> = ({ type, value, onChange, placeholder, className, ...props }) => (
  <input type={type} value={value} onChange={onChange} placeholder={placeholder} style={{ padding: '10px', margin: '5px', width: '95%' }} className={className} {...props}/>
);

interface TextAreaProps {
    value: string;
    onChange: (e: ChangeEvent<HTMLTextAreaElement>) => void;
    placeholder?: string;
}

const TextArea: FC<TextAreaProps> = ({ value, onChange, placeholder }) => (
    <textarea value={value} onChange={onChange} placeholder={placeholder} style={{ padding: '10px', margin: '5px', width: '95%', height: '100px' }} />
);

interface SelectProps {
    value: string;
    onChange: (e: ChangeEvent<HTMLSelectElement>) => void;
    children: ReactNode;
}

const Select: FC<SelectProps> = ({ value, onChange, children }) => (
  <select value={value} onChange={onChange} style={{ padding: '10px', margin: '5px', width: 'calc(95% + 14px)'}}>
    {children}
  </select>
);

interface CheckboxProps {
    checked: boolean;
    onChange: (e: ChangeEvent<HTMLInputElement>) => void;
    label: string;
}

const Checkbox: FC<CheckboxProps> = ({ checked, onChange, label }) => (
    <label style={{ margin: '5px', display: 'block' }}>
        <input type="checkbox" checked={checked} onChange={onChange} />
        {label}
    </label>
);


export default function SymphonyComparePage() {
    const globalApiKey = useApiKeyStore((state) => state.getApiKey());
    const globalProvider = useApiKeyStore((state) => state.provider);
    
    const [apiKey, setApiKey] = useState('');
    const [provider, setProvider] = useState('');
    const [tissue, setTissue] = useState('lymph node');
    const [species, setSpecies] = useState('human');
    const [celltypes, setCelltypes] = useState('B cells, T cells, Macrophages');
    const [markerSet, setMarkerSet] = useState('CD19, MS4A1, CD3D, CD8A, CD68, CD163');
    const [modelPreset, setModelPreset] = useState('symphony');
    
    // Initialize with global settings
    useEffect(() => {
        if (globalApiKey && !apiKey) {
            setApiKey(globalApiKey);
        }
        if (globalProvider && !provider) {
            setProvider(globalProvider);
        }
    }, [globalApiKey, globalProvider]);
    
    const [enableDiscussion, setEnableDiscussion] = useState(true);
    const [maxDiscussionRounds, setMaxDiscussionRounds] = useState(2);
    const [consensusThreshold, setConsensusThreshold] = useState(0.8);

    const [isLoading, setIsLoading] = useState(false);
    const [results, setResults] = useState<{reportHtml: string} | null>(null);
    const [error, setError] = useState<string | null>(null);

    const handleRunAnalysis = async () => {
        if (!apiKey || !markerSet || !celltypes) {
            setError('Please provide an API key, marker set, and cell types.');
            return;
        }
        setIsLoading(true);
        setError(null);
        setResults(null);

        try {
            const result = await symphonyCompare({
                apiKey,
                provider,
                tissue,
                species,
                celltypes: celltypes.split(',').map(s => s.trim()),
                markerSet,
                modelPreset,
                enableDiscussion,
                maxDiscussionRounds: Number(maxDiscussionRounds),
                consensusThreshold: Number(consensusThreshold),
            });
            setResults(result);
        } catch (err: any) {
            setError(`Analysis failed: ${err.message}`);
        } finally {
            setIsLoading(false);
        }
    };

    return (
        <div style={{ padding: '20px', fontFamily: 'sans-serif' }}>
            <h1>Symphony Compare</h1>
            <p>Perform multi-model cell type comparison with AI-driven consensus building.</p>

            <div style={{ display: 'grid', gridTemplateColumns: '1fr 2fr', gap: '40px', marginTop: '20px' }}>
                <div>
                    <h2>Settings</h2>
                    <Input type="text" value={apiKey} onChange={(e) => setApiKey(e.target.value)} placeholder="API Key" />
                    <Select value={provider} onChange={(e) => setProvider(e.target.value)}>
                        <option value="openrouter">OpenRouter</option>
                        <option value="openai">OpenAI</option>
                        <option value="anthropic">Anthropic</option>
                    </Select>
                    <Input type="text" value={tissue} onChange={(e) => setTissue(e.target.value)} placeholder="Tissue" />
                    <Input type="text" value={species} onChange={(e) => setSpecies(e.target.value)} placeholder="Species" />
                    
                    <TextArea value={celltypes} onChange={(e) => setCelltypes(e.target.value)} placeholder="Candidate Cell Types (comma-separated)" />
                    <TextArea value={markerSet} onChange={(e) => setMarkerSet(e.target.value)} placeholder="Marker Set (comma-separated)" />

                    <h3>Model Configuration</h3>
                     <Select value={modelPreset} onChange={(e) => setModelPreset(e.target.value)}>
                        <option value="symphony">Symphony (Claude, Gemini, Mistral)</option>
                        <option value="classic">Classic (GPT-4, GPT-3.5, Claude-2)</option>
                    </Select>

                    <h3>Discussion Settings</h3>
                    <Checkbox checked={enableDiscussion} onChange={(e) => setEnableDiscussion(e.target.checked)} label="Enable Discussion Rounds" />
                    {enableDiscussion && (
                        <>
                            <div>
                                <label>Max Discussion Rounds:</label>
                                <Input type="number" value={maxDiscussionRounds} onChange={(e) => setMaxDiscussionRounds(Number(e.target.value))} />
                            </div>
                            <div>
                                <label>Consensus Threshold: {consensusThreshold}</label>
                                <Input type="range" min="0" max="1" step="0.05" value={consensusThreshold} onChange={(e) => setConsensusThreshold(Number(e.target.value))} />
                            </div>
                        </>
                    )}

                    <Button onClick={handleRunAnalysis} disabled={isLoading}>
                        {isLoading ? 'Running Analysis...' : 'Run Symphony Compare'}
                    </Button>
                </div>

                <div>
                    <h2>Results</h2>
                    {isLoading && <p>Loading... The models are thinking.</p>}
                    {error && <p style={{ color: 'red' }}>{error}</p>}
                    {results && (
                         <div style={{ border: '1px solid #ccc', height: '80vh' }}>
                            <iframe srcDoc={results.reportHtml} style={{ width: '100%', height: '100%', border: 'none' }} title="Symphony Compare Report" />
                        </div>
                    )}
                </div>
            </div>
        </div>
    );
} 