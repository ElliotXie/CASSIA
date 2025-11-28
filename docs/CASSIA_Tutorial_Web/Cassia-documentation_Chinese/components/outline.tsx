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

interface Heading {
  id: string
  text: string
  level: number
}

export function Outline() {
  const [headings, setHeadings] = useState<Heading[]>([])
  const [activeId, setActiveId] = useState<string>("")
  const { setTheme } = useTheme()
  const githubRepoUrl = "https://github.com/ElliotXie/CASSIA"
  const biorxivUrl = "https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2"

  useEffect(() => {
    // 收集标题的函数
    const collectHeadings = () => {
      const headingElements = Array.from(document.querySelectorAll("h1, h2, h3, h4, h5, h6"))

      const headings = headingElements
        .filter((element) => element.id) // 只包含有ID的标题
        .map((element) => {
          const id = element.id
          const text = element.textContent || ""
          const level = Number.parseInt(element.tagName.substring(1))
          return { id, text, level }
        })

      setHeadings(headings)

      // 为收集的标题设置交叉观察器
      const observer = new IntersectionObserver(
        (entries) => {
          // 找到第一个正在交叉的标题
          const intersectingEntry = entries.find((entry) => entry.isIntersecting)
          if (intersectingEntry) {
            setActiveId(intersectingEntry.target.id)
          }
        },
        { rootMargin: "0px 0px -80% 0px", threshold: 0.1 },
      )

      // 观察所有标题元素
      headingElements.forEach((element) => {
        if (element.id) {
          observer.observe(element)
        }
      })

      return observer
    }

    // 初始收集
    const observer = collectHeadings()

    // 设置一个变异观察器来检测DOM的变化
    const mutationObserver = new MutationObserver(() => {
      // 清理之前的观察器
      observer.disconnect()
      // 重新收集标题
      collectHeadings()
    })

    // 开始观察文档的变化
    mutationObserver.observe(document.body, {
      childList: true,
      subtree: true,
    })

    // 清理函数
    return () => {
      observer.disconnect()
      mutationObserver.disconnect()
    }
  }, [])

  if (headings.length === 0) {
    return null
  }

  return (
    <div className="hidden lg:block">
      <div className="fixed right-8 top-20 w-56">
        <div className="mb-4 flex items-center justify-start gap-2">
          <Link href={biorxivUrl} target="_blank" rel="noopener noreferrer">
            <Button variant="outline" size="icon" aria-label="bioRxiv论文">
              <FileText className="h-[1.2rem] w-[1.2rem]" />
            </Button>
          </Link>
          <Link href={githubRepoUrl} target="_blank" rel="noopener noreferrer">
            <Button variant="outline" size="icon" aria-label="GitHub仓库">
              <Github className="h-[1.2rem] w-[1.2rem]" />
            </Button>
          </Link>
          <DropdownMenu>
            <DropdownMenuTrigger asChild>
              <Button variant="outline" size="icon">
                <Sun className="h-[1.2rem] w-[1.2rem] rotate-0 scale-100 transition-all dark:-rotate-90 dark:scale-0" />
                <Moon className="absolute h-[1.2rem] w-[1.2rem] rotate-90 scale-0 transition-all dark:rotate-0 dark:scale-100" />
                <span className="sr-only">切换主题</span>
              </Button>
            </DropdownMenuTrigger>
            <DropdownMenuContent align="end">
              <DropdownMenuItem onClick={() => setTheme("light")}>
                浅色
              </DropdownMenuItem>
              <DropdownMenuItem onClick={() => setTheme("dark")}>
                深色
              </DropdownMenuItem>
            </DropdownMenuContent>
          </DropdownMenu>
        </div>

        {headings.length > 0 && (
          <>
            <div className="text-sm font-medium text-primary/80">本页内容</div>
            <ul className="mt-2 space-y-1 text-sm">
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
          </>
        )}
      </div>
    </div>
  )
}
