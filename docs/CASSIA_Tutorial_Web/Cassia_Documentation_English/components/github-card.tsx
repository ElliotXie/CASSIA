import { Card, CardContent, CardDescription, CardFooter, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Star, GitFork, Eye } from "lucide-react"
import Link from "next/link"
import Image from "next/image"

export function GitHubCard() {
  return (
    <Card className="overflow-hidden">
      <CardHeader className="bg-slate-50 dark:bg-slate-800 pb-2">
        <div className="flex items-center gap-2">
          <Image src="/images/cassia-logo.png" alt="CASSIA Logo" width={24} height={24} className="rounded-full" />
          <CardTitle className="text-sm font-medium">ElliotXie/CASSIA</CardTitle>
        </div>
        <CardDescription className="text-xs">A multiagent LLM based single cell Annotation framework</CardDescription>
      </CardHeader>
      <CardContent className="pt-4 pb-2 text-sm">
        <p>
          CASSIA (Collaborative Agent System for Single cell Interpretable Annotation) is a tool that enhances cell type
          annotation by leveraging the power of multi-agent Large Language Models.
        </p>
      </CardContent>
      <CardFooter className="flex justify-between border-t pt-3 pb-3 bg-slate-50 dark:bg-slate-800">
        <div className="flex gap-4">
          <div className="flex items-center gap-1 text-xs">
            <Star className="h-3.5 w-3.5" />
            <span>23</span>
          </div>
          <div className="flex items-center gap-1 text-xs">
            <GitFork className="h-3.5 w-3.5" />
            <span>2</span>
          </div>
          <div className="flex items-center gap-1 text-xs">
            <Eye className="h-3.5 w-3.5" />
            <span>2</span>
          </div>
        </div>
        <Button asChild variant="outline" size="sm" className="h-7 text-xs">
          <Link href="https://github.com/ElliotXie/CASSIA" target="_blank" rel="noopener noreferrer">
            View Repository
          </Link>
        </Button>
      </CardFooter>
    </Card>
  )
}
