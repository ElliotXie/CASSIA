"use client"

import Link from "next/link"
import { usePathname } from "next/navigation"
import { useState, useEffect } from "react"
import { Menu, X, Dna } from "lucide-react"
import { cn } from "@/lib/utils"
import { Button } from "@/components/ui/button"

// 更新章节以包括分组后的CASSIA文档所有部分
const chapters = [
  {
    title: "入门指南",
    items: [
      { title: "介绍", slug: "introduction" },
      { title: "设置CASSIA", slug: "setting-up-cassia" },
    ],
  },
  {
    title: "基本CASSIA工作流程",
    items: [
      { title: "快速模式", slug: "fast-mode" },
      { title: "单簇分析", slug: "single-cluster-analysis" },
      { title: "批量处理", slug: "batch-processing" },
      { title: "质量评分和报告生成", slug: "quality-scoring-and-report-generation" },
    ],
  },
  {
    title: "高级CASSIA工作流程",
    items: [
      { title: "可选智能体简介", slug: "introduction-to-optional-agents" },
      { title: "不确定性量化（可选）", slug: "uncertainty-quantification" },
      { title: "注释增强智能体（可选）", slug: "annotation-boost" },
      { title: "比较细胞类型（可选）", slug: "compare-cell-types" },
      { title: "子聚类分析（可选）", slug: "subclustering-analysis" },
      { title: "注释增强扩展智能体（可选）", slug: "annotation-boost-extra" },
      { title: "RAG检索增强智能体（可选）", slug: "ragagent" },
    ],
  },
  {
    title: "帮助",
    items: [{ title: "故障排除", slug: "troubleshooting" }],
  },
]

export function Sidebar() {
  const pathname = usePathname()
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

  // 在移动设备上导航时关闭侧边栏
  useEffect(() => {
    if (isMobile) {
      setIsOpen(false)
    }
  }, [pathname, isMobile])

  return (
    <>
      {/* 移动设备切换按钮 */}
      <Button
        variant="ghost"
        size="icon"
        className="fixed left-4 top-4 z-50 md:hidden"
        onClick={() => setIsOpen(!isOpen)}
      >
        {isOpen ? <X className="h-5 w-5" /> : <Menu className="h-5 w-5" />}
        <span className="sr-only">切换侧边栏</span>
      </Button>

      {/* 侧边栏 */}
      <div
        className={cn(
          "fixed inset-y-0 left-0 z-40 w-64 transform bg-background/95 backdrop-blur-sm transition-transform duration-200 ease-in-out md:relative md:translate-x-0",
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
            {chapters.map((section, sectionIndex) => (
              <div key={sectionIndex} className="mb-6">
                <div className="mb-4 text-xs font-semibold uppercase tracking-wider text-muted-foreground">
                  {section.title}
                </div>
                <ul className="space-y-1">
                  {section.items.map((chapter) => (
                    <li key={chapter.slug}>
                      <Link
                        href={`/docs/${chapter.slug}`}
                        className={cn(
                          "block rounded-md px-3 py-2 text-sm font-medium transition-colors hover:bg-accent hover:text-accent-foreground",
                          pathname === `/docs/${chapter.slug}`
                            ? "bg-accent text-accent-foreground"
                            : "text-muted-foreground",
                        )}
                      >
                        {chapter.title}
                      </Link>
                    </li>
                  ))}
                </ul>
              </div>
            ))}

            <div className="mt-8 mb-4 text-xs font-semibold uppercase tracking-wider text-muted-foreground">
              资源
            </div>
            <ul className="space-y-1">
              <li>
                <Link
                  href="https://github.com/ElliotXie/CASSIA"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="block rounded-md px-3 py-2 text-sm font-medium text-muted-foreground transition-colors hover:bg-accent hover:text-accent-foreground"
                >
                  GitHub仓库
                </Link>
              </li>
              <li>
                <Link
                  href="https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="block rounded-md px-3 py-2 text-sm font-medium text-muted-foreground transition-colors hover:bg-accent hover:text-accent-foreground"
                >
                  研究论文
                </Link>
              </li>
              <li>
                <Link
                  href="/comments"
                  className="block rounded-md px-3 py-2 text-sm font-medium text-muted-foreground transition-colors hover:bg-accent hover:text-accent-foreground"
                >
                  评论区
                </Link>
              </li>
            </ul>
          </nav>
          <div className="border-t p-4">
            <div className="text-xs text-muted-foreground">CASSIA v1.2.0</div>
          </div>
        </div>
      </div>
    </>
  )
}
