"use client"

import { usePathname, useRouter } from "next/navigation"
import { useLocale, useTranslations } from "next-intl"
import { cn } from "@/lib/utils"
import { type ProgrammingLanguage, getProgrammingLanguageFromPath } from "@/contexts/programming-language-context"

interface ProgrammingLanguageSwitcherProps {
  className?: string
}

export function ProgrammingLanguageSwitcher({ className }: ProgrammingLanguageSwitcherProps) {
  const pathname = usePathname()
  const router = useRouter()
  const locale = useLocale()
  const t = useTranslations("programmingLanguage")

  // Get current programming language from URL
  const currentLang = getProgrammingLanguageFromPath(pathname)

  const switchLanguage = (newLang: ProgrammingLanguage) => {
    if (newLang === currentLang) return

    // Replace /r/ or /python/ in the pathname with the new language
    let newPathname = pathname

    if (pathname.includes("/r/")) {
      newPathname = pathname.replace("/r/", `/${newLang}/`)
    } else if (pathname.includes("/python/")) {
      newPathname = pathname.replace("/python/", `/${newLang}/`)
    } else {
      // If no language in path, we might be on a legacy URL
      // Insert language after /docs/ or /vignette/
      if (pathname.includes("/docs/")) {
        newPathname = pathname.replace("/docs/", `/docs/${newLang}/`)
      } else if (pathname.includes("/vignette/")) {
        newPathname = pathname.replace("/vignette/", `/vignette/${newLang}/`)
      }
    }

    router.push(newPathname)
  }

  return (
    <div className={cn("grid grid-cols-2 rounded-lg bg-muted p-1", className)}>
      <button
        onClick={() => switchLanguage("r")}
        className={cn(
          "px-3 py-1.5 text-sm font-medium rounded-md transition-all text-center",
          currentLang === "r"
            ? "bg-background text-foreground shadow-sm"
            : "text-muted-foreground hover:text-foreground"
        )}
      >
        {t("r")}
      </button>
      <button
        onClick={() => switchLanguage("python")}
        className={cn(
          "px-3 py-1.5 text-sm font-medium rounded-md transition-all text-center",
          currentLang === "python"
            ? "bg-background text-foreground shadow-sm"
            : "text-muted-foreground hover:text-foreground"
        )}
      >
        {t("python")}
      </button>
    </div>
  )
}
