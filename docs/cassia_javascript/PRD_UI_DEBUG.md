# CASSIA UI Debug Guide

## Overview
CASSIA web interface for cell type annotation analysis. Focus on UI debugging and component interactions.

## Key Components & Issues

### 1. Pipeline Page (`/pipeline`)
**Main Issues Fixed:**
- ✅ API key input not working → Fixed store delegation
- ✅ Provider switching not working → Fixed missing store methods  
- ✅ Example dataset loading wrong file → Updated to `batch_raw_seurat_example.csv`

**Key Files:**
- `/app/pipeline/page.tsx` - Main pipeline page
- `/components/ApiKeyInput.tsx` - API key management
- `/components/FileUpload.tsx` - File upload and example loading
- `/lib/stores/config-store.ts` - Configuration state
- `/lib/stores/api-key-store.ts` - API key state

### 2. File Upload Component
**Expected Behavior:**
- Drag & drop CSV/XLSX files
- Load example dataset button should work
- Display file metadata (rows, clusters, genes)

**Debug Points:**
- Check console for file processing errors
- Verify example file exists at `/public/examples/batch_raw_seurat_example.csv`
- Check FileUpload state updates

### 3. API Key Management
**Expected Behavior:**
- Input field should accept text
- Provider buttons (OpenRouter/Anthropic/OpenAI) should switch
- Model dropdown should update based on provider
- Values should persist in localStorage

**Debug Points:**
- Check if `useApiKeyStore` and `useConfigStore` are properly connected
- Verify localStorage persistence with DevTools → Application → Storage
- Check console for store errors

### 4. Results Display
**Expected Behavior:**
- Show real data from analysis (not mock data)
- Display actual cluster names, cell types, scores
- Real metrics: cluster count, average score, gene count

**Debug Points:**
- Verify `finalResults` structure in pipeline completion
- Check ResultsViewer CSV parsing logic
- Ensure metadata calculations are working

## Common Debug Steps

### 1. Store Issues
```javascript
// Check store state in console
useApiKeyStore.getState()
useConfigStore.getState()
useAnalysisStore.getState()
```

### 2. File Upload Issues
```javascript
// Check file processing
console.log('File processed:', data)
console.log('File metadata:', fileMetadata)
```

### 3. Pipeline Execution
```javascript
// Check pipeline progress
console.log('Pipeline state:', state)
console.log('Results:', finalResults)
```

## Development Commands
```bash
npm run dev          # Start development server
npm run build        # Test production build
npm run lint         # Check for linting errors
```

## Browser DevTools Checklist
- [ ] Console - Check for JavaScript errors
- [ ] Network - Verify API calls and file requests
- [ ] Application → Storage → Local Storage - Check persisted state
- [ ] Elements - Inspect component structure and styling

## Quick Fixes Applied
1. **API Key Store Connection**: Added delegation methods in config store
2. **Example Dataset**: Fixed file path to correct Seurat example
3. **Results Display**: Enhanced CSV parsing for real data display
4. **Dependency Conflict**: Fixed Zod version for Vercel deployment

## Test Scenarios
1. **Basic Flow**: Upload file → Enter API key → Run analysis → View results
2. **Provider Switching**: Change from OpenRouter to OpenAI, verify model updates
3. **Example Loading**: Click "Load & Download Example", verify file loads and processes
4. **Pipeline Execution**: Full end-to-end analysis with real data display