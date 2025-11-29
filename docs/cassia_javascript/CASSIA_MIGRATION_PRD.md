# CASSIA JavaScript Migration PRD

## Overview
CASSIA is being migrated from a Python package to a 100% client-side JavaScript application that runs entirely in the browser without any backend API requirements.

## Current State
- **Frontend**: Next.js application (`cassia-web/`)
- **Backend Migration**: Python functions being converted to JavaScript in `cassia-web/lib/cassia/`
- **Progress**: Some modules successfully migrated (e.g., `runCASSIA_batch.js`)

## Architecture

### Directory Structure
```
cassia-web/
├── app/                    # Next.js pages
│   ├── agents/            # Advanced analysis tools
│   │   ├── annotation-boost/   # Annotation enhancement
│   │   ├── scoring/           # Enhanced cell type scoring
│   │   ├── symphony-compare/  # Multi-model consensus analysis
│   │   └── test-symphony/     # Symphony testing interface
│   ├── batch/             # Batch processing UI
│   ├── subclustering/     # Advanced subclustering analysis
│   └── pipeline/          # Pipeline management
├── lib/
│   ├── cassia/           # Core CASSIA functions (JS & Python)
│   │   ├── runCASSIA.js         # Main analysis function
│   │   ├── runCASSIA_batch.js   # Batch processing
│   │   ├── annotationBoost.js   # Annotation enhancement
│   │   ├── scoring.js           # Cell type scoring
│   │   ├── symphonyCompare.js   # Multi-model consensus
│   │   ├── symphony_compare.py  # Python backend for Symphony
│   │   └── llm_utils.js         # LLM integration
│   └── stores/           # State management
├── components/           # React components
└── public/examples/      # Example datasets for testing
```

## Migration Status

### Completed ✅
- `runCASSIA_batch.js` - Batch cell type analysis
- Basic LLM integration (`llm_utils.js`)
- File upload/processing components
- **Symphony Compare** - Multi-model consensus analysis system
- **Enhanced Scoring Agent** - Advanced cell type scoring with batch processing
- **Subclustering Analysis** - Complete subclustering functionality with configurable parameters
- **Unified API Architecture** - Support for OpenRouter, OpenAI, and Anthropic APIs
- **Centralized Model Configuration** - Single source of truth for all model settings via `model_settings.json`
- **Smart Pipeline Presets** - Auto-configured performance/balanced modes with optimal model selection
- **Modern UI/UX** - Glass morphism design, responsive layouts, dark mode support

### In Progress 🔄
- `annotationBoost.js` - Annotation enhancement
- Performance optimizations for large datasets

### Pending ⏳
- Complete Python → JS conversion for remaining utilities
- Offline capabilities implementation

## Key Features

### 1. Symphony Compare - Multi-Model Consensus Analysis
- **Multi-model orchestration**: Leverages multiple AI models (Claude, GPT-4, Gemini) for improved accuracy
- **Consensus building**: Models engage in discussion rounds to reach agreement on cell type identification
- **Confidence scoring**: Provides reliability metrics for consensus decisions
- **Interactive HTML reports**: Downloadable analysis reports with detailed model discussions

### 2. Enhanced Cell Type Analysis
- Processes gene expression markers with advanced algorithms
- Single and batch processing modes with parallel execution
- Configurable parameters (gene count, temperature, batch size)
- Real-time progress tracking and retry mechanisms
- Quality metrics and statistical analysis

### 3. Advanced Subclustering
- Complete subclustering analysis with drag-and-drop file upload
- Configurable gene analysis range (10-200 genes)
- Temperature controls for model creativity
- Batch processing with worker thread management
- Real-time results display with table views

### 4. Data Formats & Processing
- CSV input/output with advanced parsing
- Seurat differential expression format support
- Pre-processed marker lists handling
- Multiple export formats for results
- Example datasets for testing and validation

### 5. Browser-Only Execution
- No server dependencies for core functionality
- Uses Web APIs for file handling
- Client-side CSV parsing and processing
- Hybrid architecture supporting both client-side JS and Python backend options

## Technical Details

### API Integration
- **Unified API Architecture**: Support for OpenRouter, OpenAI, and Anthropic APIs
- **Global API key management**: Centralized API key storage and reuse across tools
- **Provider abstraction**: Seamless switching between different AI service providers
- **Model configuration**: 
  - **Centralized settings**: All model configurations sourced from `model_settings.json` for consistency
  - **Provider prioritization**: OpenRouter displayed as first choice across all interfaces
  - **Smart defaults**: Auto-selects balanced preset when OpenRouter is chosen (user can override)
  - **Pipeline presets**:
    - **Performance**: Claude Sonnet 4 (annotation/boost), GPT-4o (scoring)
    - **Balanced**: Gemini 2.5 Flash (annotation/boost), Claude Haiku (scoring)
  - **User-friendly labels**: "Medium quality" instead of "low quality" for better UX
  - **Cost-tier mapping**: Models automatically categorized by performance and cost tiers
- **Robust error handling**: Automatic retry mechanisms and graceful degradation

### Parallel Processing
- **Multi-model consensus**: Orchestrates parallel execution across multiple AI models
- **Semaphore-based concurrency control**: Configurable worker thread management
- **Browser-compatible async/await**: Non-blocking processing for large datasets
- **Real-time progress tracking**: Live updates with detailed status information
- **Batch optimization**: Intelligent batching for improved performance

### Data Processing
- Marker extraction from differential expression data
- Auto-detection of column formats
- CSV generation with proper escaping

## Migration Guidelines

### When Converting Python → JavaScript:
1. Replace pandas operations with array methods
2. Use browser-compatible file handling
3. Implement proper error handling for browser environment
4. Ensure no server-side dependencies

### Testing & Validation
- **Comprehensive test suites**: Automated testing for all major features (`test-symphony/`, `test_symphony_debug.js`)
- **Debug panels**: Real-time model response inspection and validation
- **Example datasets**: Built-in example CSV files for testing and validation
- **Quality assurance workflows**: Input validation and error reporting
- **Performance monitoring**: Response time and success rate tracking
- **Manual testing UI**: Interactive testing through web interface

## Next Steps

### Phase 1: Supabase Integration & User Management
1. **User Authentication System**: Implement Supabase Auth with login/register/profile management
2. **Secure API Key Storage**: Replace localStorage with encrypted Supabase storage for API keys
3. **Analysis Results Storage**: Persistent storage for analysis history and results
4. **User Dashboard**: Personal workspace with saved analyses and settings

### Phase 2: Enhanced Data Management
1. **Analysis Sessions**: Save and resume analysis work with progress tracking
2. **Results Library**: Browse, search, and organize past analysis results
3. **Export/Import**: Enhanced data portability with cloud backup capabilities
4. **Collaboration Features**: Share results and analyses with team members

### Phase 3: Advanced Features
1. **Real-time Collaboration**: Multi-user analysis sessions with live updates
2. **Usage Analytics**: Track and analyze user patterns and model performance
3. **Advanced Visualization**: Interactive charts and plots for results
4. **API Endpoints**: REST API for programmatic access to analysis capabilities

### Phase 4: Enterprise Features
1. **Team Workspaces**: Organization-level account management
2. **Integration Plugins**: Connectors for popular bioinformatics tools
3. **Performance Monitoring**: Advanced analytics and usage tracking
4. **Offline Capabilities**: Service worker for offline processing

## Supabase Integration Architecture

### Database Schema
```sql
-- User profiles (extends Supabase Auth)
CREATE TABLE profiles (
  id UUID REFERENCES auth.users(id) PRIMARY KEY,
  email VARCHAR(255),
  full_name VARCHAR(255),
  avatar_url VARCHAR(255),
  created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
  updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Encrypted API keys storage
CREATE TABLE user_api_keys (
  id UUID DEFAULT gen_random_uuid() PRIMARY KEY,
  user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
  provider VARCHAR(50) NOT NULL, -- 'openrouter', 'anthropic', 'openai'
  encrypted_key TEXT NOT NULL,
  created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
  updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
  UNIQUE(user_id, provider)
);

-- Analysis results storage
CREATE TABLE analysis_results (
  id UUID DEFAULT gen_random_uuid() PRIMARY KEY,
  user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
  analysis_type VARCHAR(50) NOT NULL, -- 'batch', 'symphony', 'scoring', 'subclustering'
  title VARCHAR(255),
  description TEXT,
  input_data JSONB,
  results JSONB,
  settings JSONB,
  created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
  updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Analysis sessions for tracking ongoing work
CREATE TABLE analysis_sessions (
  id UUID DEFAULT gen_random_uuid() PRIMARY KEY,
  user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
  session_name VARCHAR(255),
  analysis_type VARCHAR(50),
  status VARCHAR(20) DEFAULT 'active', -- 'active', 'completed', 'archived'
  progress JSONB,
  created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
  updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);
```

### Integration Benefits
- **Secure API Key Storage**: Encrypted keys tied to user accounts instead of localStorage
- **Analysis History**: Persistent storage for all analysis results and sessions
- **Cross-Device Sync**: Access your data from anywhere with user authentication
- **Collaboration**: Share results and analyses with team members
- **Better UX**: Persistent sessions, progress tracking, and data recovery
- **Real-time Features**: Live updates and collaborative analysis capabilities

### Implementation Strategy
1. **Phase 1**: User authentication and secure API key management
2. **Phase 2**: Analysis results storage and user dashboard
3. **Phase 3**: Collaboration features and real-time updates
4. **Phase 4**: Advanced analytics and enterprise features

## Recent Achievements
✅ **Symphony Compare**: Revolutionary multi-model consensus system implemented  
✅ **Enhanced Scoring**: Advanced batch processing with quality metrics  
✅ **Subclustering**: Complete analysis workflow with configurable parameters  
✅ **Modern Architecture**: Unified API system supporting multiple providers  
✅ **Centralized Model Management**: Single source of truth via `model_settings.json` with smart defaults  
✅ **Optimized Pipeline Presets**: Performance/balanced modes with auto-selection and optimal model pairing  
✅ **Enhanced UX**: Provider prioritization, quality label improvements, and user-friendly interfaces  
✅ **Testing Framework**: Comprehensive validation and debug capabilities