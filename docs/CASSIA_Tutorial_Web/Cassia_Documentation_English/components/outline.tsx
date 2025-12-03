"use client"

import { useEffect, useState } from "react"
import { cn } from "@/lib/utils"
import { useTheme } from "next-themes"
import { Button } from "@/components/ui/button"
import { Moon, Sun, Github, FileText } from "lucide-react"
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu"
import Link from "next/link"
import { useTranslations } from "next-intl"
import { LanguageSwitcher } from "@/components/language-switcher"

interface Heading {
  id: string
  text: string
  level: number
}

export function Outline() {
  const [headings, setHeadings] = useState<Heading[]>([])
  const [activeId, setActiveId] = useState<string>("")
  const [isLargeScreen, setIsLargeScreen] = useState(false)
  const { setTheme } = useTheme()
  const t = useTranslations("outline")
  const tTheme = useTranslations("theme")
  const githubRepoUrl = "https://github.com/ElliotXie/CASSIA"
  const biorxivUrl = "https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2"

  // Check if we're on a large screen (lg breakpoint = 1024px)
  useEffect(() => {
    const checkScreenSize = () => setIsLargeScreen(window.innerWidth >= 1024)
    checkScreenSize()
    window.addEventListener('resize', checkScreenSize)
    return () => window.removeEventListener('resize', checkScreenSize)
  }, [])

  // Only run observers on large screens where outline is visible
  useEffect(() => {
    if (!isLargeScreen) return
    // Function to collect headings
    const collectHeadings = () => {
      const headingElements = Array.from(document.querySelectorAll("h1, h2, h3, h4, h5, h6"))

      const headings = headingElements
        .filter((element) => element.id) // Only include headings with IDs
        .map((element) => {
          const id = element.id
          const text = element.textContent || ""
          const level = Number.parseInt(element.tagName.substring(1))
          return { id, text, level }
        })

      setHeadings(headings)

      // Set up intersection observer for the collected headings
      const observer = new IntersectionObserver(
        (entries) => {
          // Find the first heading that's intersecting
          const intersectingEntry = entries.find((entry) => entry.isIntersecting)
          if (intersectingEntry) {
            setActiveId(intersectingEntry.target.id)
          }
        },
        { rootMargin: "0px 0px -80% 0px", threshold: 0.1 },
      )

      // Observe all heading elements
      headingElements.forEach((element) => {
        if (element.id) {
          observer.observe(element)
        }
      })

      return observer
    }

    // Initial collection
    const observer = collectHeadings()

    // Set up a mutation observer to detect changes in the DOM
    const mutationObserver = new MutationObserver(() => {
      // Clean up previous observer
      observer.disconnect()
      // Recollect headings
      collectHeadings()
    })

    // Start observing the document for changes
    mutationObserver.observe(document.body, {
      childList: true,
      subtree: true,
    })

    // Clean up function
    return () => {
      observer.disconnect()
      mutationObserver.disconnect()
    }
  }, [isLargeScreen])

  if (headings.length === 0) {
    return null
  }

  return (
    <div className="hidden lg:block">
      <div className="fixed right-8 top-20 w-56 flex flex-col max-h-[calc(100vh-100px)]">
        <div className="mb-4 flex items-center justify-start gap-2 flex-shrink-0">
          <LanguageSwitcher />
          <Link href={biorxivUrl} target="_blank" rel="noopener noreferrer">
            <Button variant="outline" size="icon" aria-label="bioRxiv Paper">
              <FileText className="h-[1.2rem] w-[1.2rem]" />
            </Button>
          </Link>
          <Link href={githubRepoUrl} target="_blank" rel="noopener noreferrer">
            <Button variant="outline" size="icon" aria-label="GitHub Repository">
              <Github className="h-[1.2rem] w-[1.2rem]" />
            </Button>
          </Link>
          <DropdownMenu>
            <DropdownMenuTrigger asChild>
              <Button variant="outline" size="icon">
                <Sun className="h-[1.2rem] w-[1.2rem] rotate-0 scale-100 transition-all dark:-rotate-90 dark:scale-0" />
                <Moon className="absolute h-[1.2rem] w-[1.2rem] rotate-90 scale-0 transition-all dark:rotate-0 dark:scale-100" />
                <span className="sr-only">{tTheme("toggleTheme")}</span>
              </Button>
            </DropdownMenuTrigger>
            <DropdownMenuContent align="end">
              <DropdownMenuItem onClick={() => setTheme("light")}>
                {tTheme("light")}
              </DropdownMenuItem>
              <DropdownMenuItem onClick={() => setTheme("dark")}>
                {tTheme("dark")}
              </DropdownMenuItem>
            </DropdownMenuContent>
          </DropdownMenu>
        </div>

        {headings.length > 0 && (
          <>
            <div className="text-sm font-medium text-primary/80 flex-shrink-0">{t("onThisPage")}</div>
            <div className="overflow-y-auto mt-2 pr-2 flex-grow pb-4">
              <ul className="space-y-1 text-sm">
          {headings.map((heading) => (
            <li key={heading.id} style={{ paddingLeft: `${(heading.level - 1) * 12}px` }}>
              <a
                href={`#${heading.id}`}
                className={cn(
                  "block py-1 text-muted-foreground transition-colors hover:text-foreground no-underline",
                  activeId === heading.id && "text-primary font-medium",
                )}
                onClick={(e) => {
                  e.preventDefault()
                  const element = document.getElementById(heading.id)
                  if (element) {
                    window.scrollTo({
                      top: element.offsetTop - 80,
                      behavior: "smooth",
                    })
                    window.history.pushState(null, "", `#${heading.id}`)
                    setActiveId(heading.id)
                    element.tabIndex = -1
                    element.focus({ preventScroll: true })
                  }
                }}
              >
                {heading.text}
              </a>
            </li>
          ))}
        </ul>
            </div>
          </>
        )}
      </div>
    </div>
  )
}
