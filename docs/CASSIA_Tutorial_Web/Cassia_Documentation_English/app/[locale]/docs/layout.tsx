'use client' // This might still be needed if other client components are used, but not for the icons anymore

import type React from "react"
import { Outline } from "@/components/outline"
// Remove imports related to icons/theme if no longer needed here
// import { useTheme } from "next-themes"
// import { Button } from "@/components/ui/button"
// import { Moon, Sun, Github, FileText } from "lucide-react"
// import {
//   DropdownMenu,
//   DropdownMenuContent,
//   DropdownMenuItem,
//   DropdownMenuTrigger,
// } from "@/components/ui/dropdown-menu"
// import Link from "next/link"

export default function DocsLayout({
  children,
}: {
  children: React.ReactNode
}) {
  // Remove theme hook and URLs if no longer needed here
  // const { setTheme } = useTheme()
  // const githubRepoUrl = "https://github.com/ElliotXie/CASSIA"
  // const biorxivUrl = "https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2"

  return (
    <div className="flex min-h-screen">
      <div className="flex-1 px-4 py-8 md:px-8 lg:px-12">
        {/* Remove Header Section */}
        {/* <div className="mb-8 flex items-center justify-end gap-2">
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
                <span className="sr-only">Toggle theme</span>
              </Button>
            </DropdownMenuTrigger>
            <DropdownMenuContent align="end">
              <DropdownMenuItem onClick={() => setTheme("light")}>
                Light
              </DropdownMenuItem>
              <DropdownMenuItem onClick={() => setTheme("dark")}>
                Dark
              </DropdownMenuItem>
            </DropdownMenuContent>
          </DropdownMenu>
        </div> */}

        {/* Main Content */}
        <div className="mx-auto max-w-3xl">{children}</div>
      </div>
      <Outline />
    </div>
  )
}
