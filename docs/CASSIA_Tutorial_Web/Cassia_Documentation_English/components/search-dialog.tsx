"use client"

import * as React from "react"
import { useRouter } from "next/navigation"
import { useLocale, useTranslations } from "next-intl"
import { FileText, Book, Search, Loader2 } from "lucide-react"
import {
  CommandDialog,
  CommandEmpty,
  CommandGroup,
  CommandInput,
  CommandItem,
  CommandList,
} from "@/components/ui/command"
import { useSearch } from "@/contexts/search-context"
import type { SearchItem } from "@/lib/search-data"

export function SearchDialog() {
  const { isOpen, setIsOpen } = useSearch()
  const router = useRouter()
  const locale = useLocale()
  const t = useTranslations("search")

  // Lazy load search data only when dialog opens
  const [searchData, setSearchData] = React.useState<SearchItem[]>([])
  const [isLoading, setIsLoading] = React.useState(false)

  React.useEffect(() => {
    if (isOpen && searchData.length === 0 && !isLoading) {
      setIsLoading(true)
      import("@/lib/search-data").then((module) => {
        setSearchData(module.searchData)
        setIsLoading(false)
      })
    }
  }, [isOpen, searchData.length, isLoading])

  const runCommand = React.useCallback((command: () => void) => {
    setIsOpen(false)
    command()
  }, [setIsOpen])

  // Filter search data by current locale and group by section
  const filteredData = searchData.filter(item => item.locale === locale || item.locale === "en")

  const docsItems = filteredData.filter(item => item.section === "docs")
  const vignetteItems = filteredData.filter(item => item.section === "vignette")

  return (
    <CommandDialog open={isOpen} onOpenChange={setIsOpen}>
      <CommandInput placeholder={t("placeholder")} />
      <CommandList>
        {isLoading ? (
          <div className="flex items-center justify-center py-6">
            <Loader2 className="h-6 w-6 animate-spin text-muted-foreground" />
          </div>
        ) : (
          <>
        <CommandEmpty>{t("noResults")}</CommandEmpty>
        <CommandGroup heading={t("docs")}>
          {docsItems.map((item) => (
            <CommandItem
              key={`${item.section}-${item.lang}-${item.slug}`}
              value={`${item.title} ${item.content}`}
              onSelect={() => {
                runCommand(() => router.push(`/${locale}${item.url.substring(3)}`))
              }}
            >
              <FileText className="mr-2 h-4 w-4" />
              <div className="flex flex-col">
                <span>{item.title}</span>
                <span className="text-xs text-muted-foreground">
                  {item.lang.toUpperCase()} - {item.content.substring(0, 60)}...
                </span>
              </div>
            </CommandItem>
          ))}
        </CommandGroup>
        <CommandGroup heading={t("vignettes")}>
          {vignetteItems.map((item) => (
            <CommandItem
              key={`${item.section}-${item.lang}-${item.slug}`}
              value={`${item.title} ${item.content}`}
              onSelect={() => {
                runCommand(() => router.push(`/${locale}${item.url.substring(3)}`))
              }}
            >
              <Book className="mr-2 h-4 w-4" />
              <div className="flex flex-col">
                <span>{item.title}</span>
                <span className="text-xs text-muted-foreground">
                  {item.lang.toUpperCase()} - {item.content.substring(0, 60)}...
                </span>
              </div>
            </CommandItem>
          ))}
        </CommandGroup>
          </>
        )}
      </CommandList>
    </CommandDialog>
  )
}

// Search button component for the sidebar
export function SearchButton() {
  const { toggle } = useSearch()
  const t = useTranslations("search")

  return (
    <button
      onClick={toggle}
      className="flex items-center gap-2 w-full px-3 py-2 text-sm text-muted-foreground rounded-md border border-input bg-background hover:bg-accent hover:text-accent-foreground transition-colors"
    >
      <Search className="h-4 w-4" />
      <span className="flex-1 text-left">{t("button")}</span>
      <kbd className="pointer-events-none hidden h-5 select-none items-center gap-1 rounded border bg-muted px-1.5 font-mono text-[10px] font-medium opacity-100 sm:flex">
        <span className="text-xs">⌘</span>K
      </kbd>
    </button>
  )
}
