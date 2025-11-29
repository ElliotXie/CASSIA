"use client";

import { useState, ChangeEvent, FC, ReactNode, useEffect } from 'react';
import { scoreAnnotationBatch } from '@/lib/cassia/scoring';
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

interface ProgressProps {
    value: number;
    max: number;
}

const Progress: FC<ProgressProps> = ({ value, max }) => (
    <div style={{ width: '100%', backgroundColor: '#f0f0f0', borderRadius: '5px', overflow: 'hidden', margin: '5px 0' }}>
        <div 
            style={{ 
                width: `${(value / max) * 100}%`, 
                backgroundColor: '#4CAF50', 
                height: '20px',
                transition: 'width 0.3s ease'
            }} 
        />
    </div>
);

export default function ScoringPage() {
    const globalApiKey = useApiKeyStore((state) => state.getApiKey());
    const globalProvider = useApiKeyStore((state) => state.provider);
    
    const [apiKey, setApiKey] = useState('');
    const [provider, setProvider] = useState('');
    const [model, setModel] = useState('deepseek/deepseek-chat-v3-0324');
    const [maxWorkers, setMaxWorkers] = useState(4);
    const [maxRetries, setMaxRetries] = useState(1);
    
    // Initialize with global settings
    useEffect(() => {
        if (globalApiKey && !apiKey) {
            setApiKey(globalApiKey);
        }
        if (globalProvider && !provider) {
            setProvider(globalProvider);
        }
    }, [globalApiKey, globalProvider]);
    
    const [csvFile, setCsvFile] = useState<File | null>(null);
    const [csvData, setCsvData] = useState<string>('');
    
    const [isLoading, setIsLoading] = useState(false);
    const [progress, setProgress] = useState({ completed: 0, total: 0, percentage: 0 });
    const [results, setResults] = useState<any>(null);
    const [error, setError] = useState<string | null>(null);

    const handleFileUpload = (e: ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (file) {
            setCsvFile(file);
            setError(null);
            const reader = new FileReader();
            reader.onload = (event) => {
                const text = event.target?.result as string;
                setCsvData(text);
                
                // Validate CSV format
                try {
                    const lines = text.trim().split('\n');
                    if (lines.length < 2) {
                        setError('CSV file appears to be empty or has no data rows.');
                        return;
                    }
                    
                    const headers = lines[0].split(',').map(h => h.trim().replace(/"/g, ''));
                    
                    // Check for conversation history column
                    const hasConversationHistory = headers.some(h => 
                        h.toLowerCase().includes('conversation') && h.toLowerCase().includes('history')
                    );
                    
                    if (!hasConversationHistory) {
                        setError(`This appears to be a summary CSV without conversation history. Please upload the "*_full.csv" file from runCASSIA batch analysis. Found columns: ${headers.join(', ')}`);
                        return;
                    }
                    
                    // Check for required columns
                    const requiredColumns = ['species', 'tissue'];
                    const missingColumns = requiredColumns.filter(req => 
                        !headers.some(h => h.toLowerCase().includes(req))
                    );
                    
                    if (missingColumns.length > 0) {
                        setError(`Missing required columns: ${missingColumns.join(', ')}. This may not be a runCASSIA batch output file.`);
                        return;
                    }
                    
                    console.log(`✅ CSV validation passed. Found ${lines.length - 1} data rows with conversation history.`);
                } catch (validationError) {
                    console.warn('CSV validation failed:', validationError);
                }
            };
            reader.readAsText(file);
        }
    };

    const handleRunScoring = async () => {
        if (!apiKey) {
            setError('Please provide an API key.');
            return;
        }
        
        if (!csvData) {
            setError('Please upload a CSV file with annotation results.');
            return;
        }
        
        setIsLoading(true);
        setError(null);
        setResults(null);
        setProgress({ completed: 0, total: 0, percentage: 0 });

        try {
            const result = await scoreAnnotationBatch({
                csvData,
                apiKey,
                maxWorkers: Number(maxWorkers),
                model,
                provider,
                maxRetries: Number(maxRetries),
                onProgress: (progressData: any) => {
                    setProgress(progressData);
                }
            } as any);
            
            setResults(result);
            console.log('Scoring completed:', result);
        } catch (err: any) {
            setError(`Scoring failed: ${err.message}`);
        } finally {
            setIsLoading(false);
        }
    };

    const downloadResults = () => {
        if (!results?.csvContent) return;
        
        const blob = new Blob([results.csvContent], { type: 'text/csv' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `${csvFile?.name.replace('.csv', '') || 'results'}_scored.csv`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    };

    const loadExampleData = async () => {
        try {
            setError(null);
            const response = await fetch('/examples/scoring_example_full.csv');
            const text = await response.text();
            setCsvData(text);
            
            // Create a mock file object for display
            const mockFile = new File([text], 'scoring_example_full.csv', { type: 'text/csv' });
            setCsvFile(mockFile);
            
            console.log('✅ Example data loaded successfully');
        } catch (error) {
            setError('Failed to load example data. Please try uploading your own file.');
        }
    };

    return (
        <div style={{ padding: '20px', fontFamily: 'sans-serif' }}>
            <h1>CASSIA Scoring Agent</h1>
            <p>Score cell type annotation results using AI evaluation. Upload the <strong>full CSV output from runCASSIA batch analysis</strong> (containing the "Conversation History" column) to get quality scores and reasoning for each annotation.</p>

            <div style={{ display: 'grid', gridTemplateColumns: '1fr 2fr', gap: '40px', marginTop: '20px' }}>
                <div>
                    <h2>Settings</h2>
                    
                    <div>
                        <label>API Key:</label>
                        <Input 
                            type="password" 
                            value={apiKey} 
                            onChange={(e) => setApiKey(e.target.value)} 
                            placeholder="Enter your API key" 
                        />
                    </div>
                    
                    <div>
                        <label>Provider:</label>
                        <Select value={provider} onChange={(e) => setProvider(e.target.value)}>
                            <option value="openrouter">OpenRouter</option>
                            <option value="openai">OpenAI</option>
                            <option value="anthropic">Anthropic</option>
                        </Select>
                    </div>
                    
                    <div>
                        <label>Model:</label>
                        <Input 
                            type="text" 
                            value={model} 
                            onChange={(e) => setModel(e.target.value)} 
                            placeholder="Model name" 
                        />
                    </div>
                    
                    <div>
                        <label>Max Workers: {maxWorkers}</label>
                        <Input 
                            type="range" 
                            min="1" 
                            max="10" 
                            step="1" 
                            value={maxWorkers} 
                            onChange={(e) => setMaxWorkers(Number(e.target.value))} 
                        />
                    </div>
                    
                    <div>
                        <label>Max Retries: {maxRetries}</label>
                        <Input 
                            type="range" 
                            min="0" 
                            max="5" 
                            step="1" 
                            value={maxRetries} 
                            onChange={(e) => setMaxRetries(Number(e.target.value))} 
                        />
                    </div>

                    <h2>Data Input</h2>
                    
                    <div>
                        <label>Upload CSV File:</label>
                        <input 
                            type="file" 
                            accept=".csv" 
                            onChange={handleFileUpload}
                            className="w-full px-3 py-2 border rounded-md"
                        />
                        <div style={{ marginTop: '5px' }}>
                            <button 
                                onClick={loadExampleData} 
                                disabled={false}
                                className="px-2 py-1 text-xs bg-blue-600 text-white rounded hover:bg-blue-700"
                            >
                                Load Example Data
                            </button>
                        </div>
                        {csvFile && <p style={{ fontSize: '12px', color: '#666' }}>
                            File: {csvFile.name} ({(csvFile.size / 1024).toFixed(1)}KB)
                        </p>}
                    </div>
                    
                    <div style={{ marginTop: '20px' }}>
                        <div style={{ fontSize: '14px', color: '#666', backgroundColor: '#f9f9f9', padding: '10px', borderRadius: '5px' }}>
                            <strong>Expected Input:</strong> Full CSV output from runCASSIA batch analysis
                            <br />
                            <strong>Required columns:</strong>
                            <ul style={{ margin: '5px 0', paddingLeft: '20px' }}>
                                <li>Species</li>
                                <li>Tissue</li>
                                <li>Marker List (gene markers)</li>
                                <li><strong>Conversation History</strong> (the annotation conversation to be scored)</li>
                            </ul>
                            <em>Note: Use the "*_full.csv" file from batch analysis, not the summary version.</em>
                        </div>
                    </div>

                    <Button onClick={handleRunScoring} disabled={isLoading || !csvData || !apiKey}>
                        {isLoading ? 'Running Scoring...' : 'Start Scoring'}
                    </Button>
                </div>

                <div>
                    <h2>Results</h2>
                    
                    {isLoading && (
                        <div>
                            <p>Scoring in progress... {progress.percentage}% complete</p>
                            <Progress value={progress.completed} max={progress.total} />
                            <p style={{ fontSize: '12px', color: '#666' }}>
                                {progress.completed} of {progress.total} rows processed
                            </p>
                        </div>
                    )}
                    
                    {error && (
                        <div style={{ color: 'red', padding: '10px', backgroundColor: '#ffe6e6', borderRadius: '5px' }}>
                            {error}
                        </div>
                    )}
                    
                    {results && (
                        <div style={{ border: '1px solid #ccc', padding: '15px', borderRadius: '5px' }}>
                            <h3>Scoring Complete!</h3>
                            <div style={{ marginBottom: '15px' }}>
                                <p><strong>Total Rows:</strong> {results.totalRows}</p>
                                <p><strong>Processed Rows:</strong> {results.processedRows}</p>
                                <p><strong>Already Scored:</strong> {results.alreadyScored}</p>
                            </div>
                            
                            <Button onClick={downloadResults} disabled={false}>
                                Download Scored Results
                            </Button>
                            
                            <div style={{ marginTop: '20px' }}>
                                <h4>Sample Results:</h4>
                                <div style={{ maxHeight: '400px', overflow: 'auto', fontSize: '12px' }}>
                                    {results.results.slice(0, 5).map((row: any, idx: number) => (
                                        <div key={idx} style={{ border: '1px solid #eee', margin: '5px 0', padding: '10px' }}>
                                            <p><strong>Species:</strong> {row.Species}</p>
                                            <p><strong>Tissue:</strong> {row.Tissue}</p>
                                            <p><strong>Score:</strong> {row.Score}</p>
                                            <p><strong>Reasoning:</strong> {String(row.Scoring_Reasoning).substring(0, 200)}...</p>
                                        </div>
                                    ))}
                                    {results.results.length > 5 && (
                                        <p style={{ fontStyle: 'italic' }}>
                                            ... and {results.results.length - 5} more rows
                                        </p>
                                    )}
                                </div>
                            </div>
                        </div>
                    )}
                </div>
            </div>
        </div>
    );
}