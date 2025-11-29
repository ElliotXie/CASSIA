# Product Requirements Document (PRD): CASSIA Web Application

## Overview

**Product Name**: CASSIA (Cell Annotation using Single-cell Sequencing Intelligence Assistant)  
**Version**: 1.0.0  
**Date**: January 2025  
**Status**: In Development - Data Persistence Phase

## Executive Summary

CASSIA is an AI-powered web application that automates cell type annotation for single-cell RNA-seq analysis. It leverages advanced language models to transform biological marker data into comprehensive cell type annotations, achieving 15-40% better accuracy than traditional methods while processing analyses in under 2 minutes.

## Current Architecture

### Frontend Stack
- **Framework**: Next.js 15.3.5 with TypeScript
- **UI Library**: Tailwind CSS with custom glassmorphism design
- **State Management**: Zustand with persistence
- **Form Handling**: React Hook Form with Zod validation
- **File Processing**: PapaParse (CSV), XLSX (Excel)

### Backend Services
- **AI Providers**: OpenRouter, OpenAI, Anthropic
- **Authentication**: Supabase Auth (fully integrated)
- **Database**: Supabase PostgreSQL (schema deployed, API key storage active)
- **Contact System**: EmailJS

### Application Structure

```
cassia-web/
├── app/                    # Next.js app router pages
│   ├── agents/            # Specialized analysis tools
│   │   ├── annotation-boost/
│   │   ├── merging/
│   │   ├── scoring/
│   │   └── symphony-compare/
│   ├── batch/             # Batch processing
│   ├── pipeline/          # Full pipeline
│   └── subclustering/     # Subclustering analysis
├── components/            # React components
│   ├── auth/             # Authentication UI (fully functional)
│   ├── dashboard/        # User dashboard (connected to Supabase)
│   └── ui/              # Reusable UI components
├── lib/                  # Core libraries
│   ├── cassia/          # Analysis algorithms
│   ├── stores/          # Zustand stores
│   └── supabase/        # Supabase configuration
└── public/              # Static assets and examples
```

## Core Features

### 1. Analysis Modules

#### CASSIA Pipeline (Main Feature)
- Full automated annotation workflow
- Quality scoring and validation
- Annotation merging for consensus
- Boost analysis for refinement
- HTML report generation

#### CASSIA Batch
- Parallel processing for multiple clusters
- Configurable worker threads
- Progress tracking
- Bulk result download

#### Specialized Agents
- **Annotation Boost**: Enhances initial annotations
- **Symphony Compare**: Multi-model comparison
- **Scoring Agent**: Quality assessment
- **Merging Agent**: Consensus building
- **Subclustering**: Hierarchical analysis

### 2. User Interface

- Modern, responsive design with glassmorphism effects
- 4-step quick start guide
- Drag-and-drop file upload
- Real-time progress tracking
- Interactive results viewer
- API key configuration panel

## Authentication Model

### Open Access Design

CASIA follows an **open access model** where all core analysis features are available to both authenticated and non-authenticated users. This design ensures maximum accessibility while providing enhanced features for registered users.

### Non-Authenticated Users

**Full Access To:**
- ✅ All analysis modules (Pipeline, Batch, Agents)
- ✅ File upload and processing
- ✅ Real-time results viewing
- ✅ Analysis report generation
- ✅ API key configuration (stored locally)
- ✅ All example files and documentation

**Limitations:**
- ❌ API keys stored only in browser localStorage
- ❌ No analysis history persistence
- ❌ No cross-device access to saved settings
- ❌ Results lost on browser refresh/clear

### Authenticated Users

**All Non-Authenticated Features Plus:**
- ✅ **Load API Keys Button** - Retrieve saved API keys from account
- ✅ Encrypted API key storage in Supabase
- ✅ Cross-device API key synchronization
- ❌ Persistent analysis history (not yet implemented)
- ❌ User dashboard with saved results (UI ready, data not connected)
- ❌ Profile management (UI ready, functionality pending)
- ❌ Analysis session tracking (not yet implemented)

## Current State vs. Target State

### Authentication & User Management

**Current State**:
- ✅ Open access model implemented
- ✅ Full functionality for non-authenticated users
- ✅ Auth UI components built and functional
- ✅ Supabase Auth fully configured and working
- ✅ Load API Keys button for authenticated users
- ✅ Sign-up/sign-in flows working
- ✅ Session persistence working
- ✅ Professional error notification modal for authentication errors

**Target State**:
- ✅ Continue open access model
- ❌ Enhanced user profiles with preferences
- ✅ Secure session management (basic implementation done)
- ❌ Password reset flow (not yet implemented)

### Data Storage

**Current State**:
- ✅ Dual API key storage (localStorage + Supabase) - fully working
- ✅ Load API Keys button for authenticated users - fully working
- ✅ API key encryption and decryption - implemented
- ❌ Results not persisted automatically
- ❌ Limited analysis history
- ✅ Database schema deployed and active

**Target State**:
- ✅ Encrypted API key storage in Supabase (completed)
- ❌ Persistent analysis results (in progress)
- ❌ Analysis history and retrieval (planned)
- ✅ Cross-device API key synchronization (completed)

### Database Schema (Ready for Integration)

```sql
-- Core tables defined in supabase_schema.sql
├── profiles (user profiles)
├── api_keys (encrypted storage)
├── analysis_results (job history)
├── analysis_sessions (usage tracking)
└── analysis_files (file metadata)
```

## Supabase Integration Roadmap

### Phase 1: Authentication (Completed)
- [x] Maintain open access model
- [x] Implement sign-up/sign-in flows
- [x] Add session persistence
- [x] Create Load API Keys functionality
- [x] Implement API key encryption/decryption
- [x] Deploy database schema

### Phase 2: Data Persistence (Current Phase)
- [x] Dual API key storage (localStorage + Supabase)
- [ ] Save analysis results to database
- [ ] Implement file upload to Supabase Storage
- [ ] Add analysis history retrieval
- [ ] Connect dashboard to analysis data

### Phase 3: User Features (Planned)
- [ ] Connect user dashboard to real analysis data
- [ ] Implement profile management functionality
- [ ] Add usage analytics
- [ ] Create API key management UI (beyond Load button)
- [ ] Implement password reset flow

### Phase 4: Advanced Features
- [ ] Team collaboration
- [ ] Result sharing
- [ ] Usage quotas
- [ ] Billing integration

## Technical Considerations

### Security
- Row Level Security (RLS) policies implemented
- API keys encrypted at rest
- Secure session management
- CORS configuration for API calls

### Performance
- Optimized database indexes
- Client-side caching with Zustand
- Lazy loading for large results
- Background job processing

### Scalability
- Stateless architecture
- Database connection pooling
- CDN for static assets
- Horizontal scaling ready

## Success Metrics

- User adoption rate
- Analysis completion rate
- Average processing time < 2 minutes
- User retention (30-day)
- API reliability > 99.9%

## Future Enhancements

1. **Advanced Analytics**
   - Comparative analysis across datasets
   - Longitudinal study support
   - Custom marker libraries

2. **Collaboration Features**
   - Team workspaces
   - Shared analysis templates
   - Commenting system

3. **Integration Ecosystem**
   - REST API for programmatic access
   - Python/R SDK
   - Jupyter notebook extension

## Development Notes for Future Agents

### Quick Start for Supabase Integration

1. **Environment Setup**
   ```bash
   # Required environment variables
   NEXT_PUBLIC_SUPABASE_URL=your_supabase_url
   NEXT_PUBLIC_SUPABASE_ANON_KEY=your_anon_key
   ```

2. **Key Files to Review**
   - `lib/supabase/client.ts` - Supabase client setup
   - `components/auth/AuthProvider.tsx` - Auth context
   - `lib/stores/auth-store.ts` - Auth state management
   - `lib/stores/api-key-store-simple.ts` - API key storage implementation
   - `components/LoadApiKeysButton.tsx` - Load API keys functionality
   - `supabase_schema.sql` - Database structure

3. **Integration Priority**
   - Maintain open access model
   - Enhance authenticated user features
   - Focus on result persistence and history

### Current Status

1. **Open Access Working** - All routes are public by design ✅
2. **Authentication System** - Sign-up/sign-in flows fully functional ✅
3. **Dual API Key Storage** - localStorage for all users, Supabase for authenticated ✅
4. **Load API Keys Feature** - Working button to retrieve saved keys ✅
5. **Database Schema** - Deployed and active for API key storage ✅
6. **Dashboard UI** - Built but not connected to analysis data ⚠️

### Next Priority Tasks

1. **Analysis Result Persistence** - Save analysis results to database
2. **Connect Dashboard Data** - Display real analysis history in dashboard
3. **Profile Management** - Enable users to update profile information
4. **Password Reset Flow** - Implement forgot password functionality
5. **Usage Analytics** - Track user engagement and analysis patterns

### Testing Considerations

- ✅ Auth flows working (sign-up/sign-in tested)
- ✅ Dual storage system (localStorage + Supabase tested)
- ✅ Load API Keys functionality tested
- ❌ Performance testing with Supabase (needs comprehensive testing)
- ❌ Security audit for RLS policies (needs review)

## Contact

For questions or support regarding CASSIA development:
- GitHub Issues: [Repository URL]
- Email: [Contact Email]
- Documentation: See README.md and SUPABASE_SETUP.md