"use client"

import { createContext, useContext, useState, useEffect, ReactNode } from "react"

export type ProgrammingLanguage = "r" | "python"

interface ProgrammingLanguageContextType {
  language: ProgrammingLanguage
  setLanguage: (lang: ProgrammingLanguage) => void
}

const ProgrammingLanguageContext = createContext<ProgrammingLanguageContextType | undefined>(undefined)

const STORAGE_KEY = "cassia-programming-language"

export function ProgrammingLanguageProvider({ children }: { children: ReactNode }) {
  const [language, setLanguageState] = useState<ProgrammingLanguage>("r")
  const [isHydrated, setIsHydrated] = useState(false)

  // Load preference from localStorage on mount
  useEffect(() => {
    const stored = localStorage.getItem(STORAGE_KEY)
    if (stored === "r" || stored === "python") {
      setLanguageState(stored)
    }
    setIsHydrated(true)
  }, [])

  // Save preference to localStorage when it changes
  const setLanguage = (lang: ProgrammingLanguage) => {
    setLanguageState(lang)
    localStorage.setItem(STORAGE_KEY, lang)
  }

  // Prevent hydration mismatch by not rendering until client-side
  if (!isHydrated) {
    return null
  }

  return (
    <ProgrammingLanguageContext.Provider value={{ language, setLanguage }}>
      {children}
    </ProgrammingLanguageContext.Provider>
  )
}

export function useProgrammingLanguage() {
  const context = useContext(ProgrammingLanguageContext)
  if (context === undefined) {
    throw new Error("useProgrammingLanguage must be used within a ProgrammingLanguageProvider")
  }
  return context
}

// Helper to extract programming language from URL path
export function getProgrammingLanguageFromPath(pathname: string): ProgrammingLanguage {
  if (pathname.includes("/python/")) {
    return "python"
  }
  return "r"
}
