"use client"

import { usePathname } from "next/navigation"
import { useState, useEffect } from "react"
import { Menu, X, Dna } from "lucide-react"
import { cn } from "@/lib/utils"
import { Button } from "@/components/ui/button"
import { useTranslations, useLocale } from "next-intl"
import { Link } from "@/i18n/routing"

// Define chapter structure with translation keys
const chapters = [
  {
    titleKey: "gettingStarted",
    items: [
      { titleKey: "introduction", slug: "introduction" },
      { titleKey: "settingUpCassia", slug: "setting-up-cassia" },
    ],
  },
  {
    titleKey: "basicWorkflow",
    items: [
      { titleKey: "fastMode", slug: "fast-mode" },
      { titleKey: "singleClusterAnalysis", slug: "single-cluster-analysis" },
      { titleKey: "batchProcessing", slug: "batch-processing" },
      { titleKey: "qualityScoring", slug: "quality-scoring-and-report-generation" },
    ],
  },
  {
    titleKey: "advancedWorkflow",
    items: [
      { titleKey: "introToOptionalAgents", slug: "introduction-to-optional-agents" },
      { titleKey: "uncertaintyQuantification", slug: "uncertainty-quantification" },
      { titleKey: "annotationBoost", slug: "annotation-boost" },
      { titleKey: "compareCellTypes", slug: "compare-cell-types" },
      { titleKey: "subclusteringAnalysis", slug: "subclustering-analysis" },
      { titleKey: "annotationBoostExtra", slug: "annotation-boost-extra" },
      { titleKey: "ragAgent", slug: "ragagent" },
    ],
  },
  {
    titleKey: "help",
    items: [{ titleKey: "troubleshooting", slug: "troubleshooting" }],
  },
]

// Define vignette chapters
const vignetteChapters = [
  {
    titleKey: "gettingStarted",
    items: [
      { titleKey: "introduction", slug: "introduction" },
    ],
  },
  {
    titleKey: "analysisWithMarker",
    items: [
      { titleKey: "basicAnnotation", slug: "basic-annotation" },
      { titleKey: "extendedAnalysis", slug: "extended-analysis" },
    ],
  },
  {
    titleKey: "analysisWithSeurat",
    items: [
      { titleKey: "clusteringAndAnnotation", slug: "clustering-and-annotation" },
      { titleKey: "fullWorkflow", slug: "full-workflow-best-practices" }
    ],
  }
]


export function Sidebar() {
  const pathname = usePathname()
  const locale = useLocale()
  const tNav = useTranslations("navigation")
  const tSidebar = useTranslations("sidebar")
  const tDocs = useTranslations("docs")
  const tVignettes = useTranslations("vignettes")
  const tCommon = useTranslations("common")

  const [isOpen, setIsOpen] = useState(false)
  const [isMobile, setIsMobile] = useState(false)

  useEffect(() => {
    const checkMobile = () => {
      setIsMobile(window.innerWidth < 768)
    }

    checkMobile()
    window.addEventListener("resize", checkMobile)

    return () => {
      window.removeEventListener("resize", checkMobile)
    }
  }, [])

  useEffect(() => {
    if (isMobile) {
      setIsOpen(false)
    } else {
      setIsOpen(true)
    }
  }, [isMobile])

  // Close sidebar when navigating on mobile
  useEffect(() => {
    if (isMobile) {
      setIsOpen(false)
    }
  }, [pathname, isMobile])

  // Determine if we're in the vignette section (accounting for locale prefix)
  const isVignettePath = pathname?.includes("/vignette")

  // Select which sections to show based on the current path
  const sectionsToShow = isVignettePath ? vignetteChapters : chapters

  // Get the section title translation
  const getSectionTitle = (titleKey: string) => {
    return tSidebar(titleKey as any)
  }

  // Get the item title translation
  const getItemTitle = (titleKey: string, isVignette: boolean) => {
    if (isVignette) {
      return tVignettes(titleKey as any)
    }
    return tDocs(titleKey as any)
  }

  return (
    <>
      {/* Mobile toggle button */}
      <Button
        variant="ghost"
        size="icon"
        className="fixed left-4 top-4 z-50 md:hidden"
        onClick={() => setIsOpen(!isOpen)}
      >
        {isOpen ? <X className="h-5 w-5" /> : <Menu className="h-5 w-5" />}
        <span className="sr-only">Toggle sidebar</span>
      </Button>

      {/* Sidebar */}
      <div
        className={cn(
          "fixed inset-y-0 left-0 z-40 w-64 transform bg-background/95 backdrop-blur-sm transition-transform duration-200 ease-in-out md:sticky md:top-0 md:h-screen md:translate-x-0",
          isOpen ? "translate-x-0" : "-translate-x-full",
        )}
      >
        <div className="flex h-full flex-col border-r">
          <div className="flex h-20 items-center border-b px-6">
            <Link href="/" className="flex items-center gap-2">
              <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-primary text-primary-foreground">
                <Dna className="h-6 w-6" />
              </div>
              <span className="font-semibold text-xl tracking-tight">CASSIA</span>
            </Link>
          </div>
          <nav className="flex-1 overflow-y-auto p-6">
            {/* Navigation tabs with equal width and centered content */}
            <div className="mb-6 grid grid-cols-2 border-b">
              <div className="flex justify-center">
                <Link
                  href="/docs/introduction"
                  className={cn(
                    "px-4 py-2 text-sm font-medium flex-1 text-center",
                    !isVignettePath ? "border-b-2 border-primary text-primary" : "text-muted-foreground"
                  )}
                >
                  {tNav("docs")}
                </Link>
              </div>
              <div className="flex justify-center">
                <Link
                  href="/vignette/introduction"
                  className={cn(
                    "px-4 py-2 text-sm font-medium flex-1 text-center",
                    isVignettePath ? "border-b-2 border-primary text-primary" : "text-muted-foreground"
                  )}
                >
                  {tNav("vignettes")}
                </Link>
              </div>
            </div>

            {sectionsToShow.map((section, sectionIndex) => (
              <div key={sectionIndex} className="mb-6">
                <div className="mb-4 text-xs font-semibold uppercase tracking-wider text-muted-foreground">
                  {getSectionTitle(section.titleKey)}
                </div>
                <ul className="space-y-1">
                  {section.items.map((chapter) => (
                    <li key={chapter.slug}>
                      <Link
                        href={`${isVignettePath ? "/vignette" : "/docs"}/${chapter.slug}`}
                        className={cn(
                          "block rounded-md px-3 py-2 text-sm font-medium transition-colors hover:bg-accent hover:text-accent-foreground",
                          pathname?.endsWith(`${isVignettePath ? "/vignette" : "/docs"}/${chapter.slug}`)
                            ? "bg-accent text-accent-foreground"
                            : "text-muted-foreground",
                        )}
                      >
                        {getItemTitle(chapter.titleKey, isVignettePath)}
                      </Link>
                    </li>
                  ))}
                </ul>
              </div>
            ))}

            <div className="mt-8 mb-4 text-xs font-semibold uppercase tracking-wider text-muted-foreground">
              {tNav("resources")}
            </div>
            <ul className="space-y-1">
              <li>
                <a
                  href="https://github.com/ElliotXie/CASSIA"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="block rounded-md px-3 py-2 text-sm font-medium text-muted-foreground transition-colors hover:bg-accent hover:text-accent-foreground"
                >
                  {tNav("githubRepository")}
                </a>
              </li>
              <li>
                <a
                  href="https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="block rounded-md px-3 py-2 text-sm font-medium text-muted-foreground transition-colors hover:bg-accent hover:text-accent-foreground"
                >
                  {tNav("researchPaper")}
                </a>
              </li>
              <li>
                <Link
                  href="/comments"
                  className="block rounded-md px-3 py-2 text-sm font-medium text-muted-foreground transition-colors hover:bg-accent hover:text-accent-foreground"
                >
                  {tNav("commentSection")}
                </Link>
              </li>
            </ul>
          </nav>
          <div className="border-t p-4">
            <div className="text-xs text-muted-foreground">{tCommon("version")}</div>
          </div>
        </div>
      </div>
    </>
  )
}
