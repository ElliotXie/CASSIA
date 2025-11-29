/**
 * Utility to clean up old localStorage entries that might cause migration issues
 */

export const cleanupOldStorage = () => {
  if (typeof window === 'undefined') return
  
  try {
    // List of old storage keys that might cause conflicts
    const oldKeys = [
      'cassia-api-key-storage',
      'cassia-auth-storage',
      'cassia-results-storage'
    ]
    
    let cleaned = false
    
    oldKeys.forEach(key => {
      const stored = localStorage.getItem(key)
      if (stored) {
        try {
          const parsed = JSON.parse(stored)
          // Check if this is an old version that might cause issues
          if (parsed.version === 1 && key === 'cassia-api-key-storage') {
            console.log(`Cleaning up old storage format for ${key}`)
            localStorage.removeItem(key)
            cleaned = true
          }
        } catch (e) {
          // If we can't parse it, it's probably corrupted, so remove it
          console.log(`Removing corrupted storage for ${key}`)
          localStorage.removeItem(key)
          cleaned = true
        }
      }
    })
    
    if (cleaned) {
      console.log('Storage cleanup completed. Please refresh the page.')
      return true
    }
    
    return false
  } catch (error) {
    console.error('Error during storage cleanup:', error)
    return false
  }
}

/**
 * Force clean all CASSIA-related storage (use with caution)
 */
export const forceCleanStorage = () => {
  if (typeof window === 'undefined') return
  
  try {
    const keys = Object.keys(localStorage)
    const cassiaKeys = keys.filter(key => key.startsWith('cassia-'))
    
    cassiaKeys.forEach(key => {
      localStorage.removeItem(key)
      console.log(`Removed ${key} from localStorage`)
    })
    
    console.log('Force cleanup completed. All CASSIA storage cleared.')
  } catch (error) {
    console.error('Error during force cleanup:', error)
  }
}

/**
 * Check storage health and report issues
 */
export const checkStorageHealth = () => {
  if (typeof window === 'undefined') return
  
  const report = {
    hasOldFormat: false,
    hasNewFormat: false,
    hasConflicts: false,
    keys: [] as string[]
  }
  
  try {
    const keys = Object.keys(localStorage)
    const cassiaKeys = keys.filter(key => key.startsWith('cassia-'))
    
    report.keys = cassiaKeys
    
    cassiaKeys.forEach(key => {
      const stored = localStorage.getItem(key)
      if (stored) {
        try {
          const parsed = JSON.parse(stored)
          if (parsed.version === 1) {
            report.hasOldFormat = true
          } else if (parsed.version === 2) {
            report.hasNewFormat = true
          }
        } catch (e) {
          // Corrupted storage
          report.hasConflicts = true
        }
      }
    })
    
    if (report.hasOldFormat && report.hasNewFormat) {
      report.hasConflicts = true
    }
    
    console.log('Storage Health Report:', report)
    return report
  } catch (error) {
    console.error('Error checking storage health:', error)
    return report
  }
}