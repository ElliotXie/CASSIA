import Link from "next/link"
import { Button } from "@/components/ui/button"
import { ArrowRight, Dna, Microscope, FileSearch, CheckCircle } from "lucide-react"
import { GitHubCard } from "@/components/github-card"

export default function Home() {
  return (
    <div className="flex min-h-screen flex-col">
      <div className="flex flex-col items-center justify-center px-4 py-20 text-center md:py-32">
        <div className="flex h-20 w-20 items-center justify-center rounded-2xl bg-primary text-primary-foreground mb-6">
          <Dna className="h-12 w-12" />
        </div>
        <h1 className="text-4xl font-extrabold tracking-tight sm:text-5xl md:text-6xl lg:text-7xl font-mono text-transparent bg-clip-text bg-gradient-to-r from-blue-500 to-cyan-500">CASSIA</h1>
        <p className="mt-4 text-xl text-muted-foreground max-w-2xl">
          为研究人员和生物信息学家提供的高级单细胞分析平台
        </p>
        <div className="mt-8 flex flex-wrap justify-center gap-4">
          <Button asChild size="lg" className="h-12 px-8">
            <Link href="/docs/introduction">
              开始使用
              <ArrowRight className="ml-2 h-4 w-4" />
            </Link>
          </Button>
          <Button asChild variant="outline" size="lg" className="h-12 px-8">
            <Link href="https://github.com/ElliotXie/CASSIA" target="_blank" rel="noopener noreferrer">
              GitHub仓库
            </Link>
          </Button>
          <Button asChild variant="outline" size="lg" className="h-12 px-8">
            <Link
              href="https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2"
              target="_blank"
              rel="noopener noreferrer"
            >
              阅读论文
            </Link>
          </Button>
        </div>
      </div>

      <div className="container mx-auto px-4 py-12 md:py-24">
        <div className="grid gap-8 md:grid-cols-3">
          <div className="rounded-lg border bg-card p-6 shadow-sm transition-all hover:shadow-md">
            <div className="mb-4 flex h-10 w-10 items-center justify-center rounded-full bg-primary/10 text-primary">
              <Microscope className="h-5 w-5" />
            </div>
            <h3 className="text-lg font-semibold">可解释性</h3>
            <p className="mt-2 text-muted-foreground">
              为每个注释提供详细的推理过程，使分析过程透明且易于理解。
            </p>
          </div>
          <div className="rounded-lg border bg-card p-6 shadow-sm transition-all hover:shadow-md">
            <div className="mb-4 flex h-10 w-10 items-center justify-center rounded-full bg-primary/10 text-primary">
              <FileSearch className="h-5 w-5" />
            </div>
            <h3 className="text-lg font-semibold">无参考注释</h3>
            <p className="mt-2 text-muted-foreground">
              无需参考数据集即可注释细胞类型，能够发现新的细胞群体。
            </p>
          </div>
          <div className="rounded-lg border bg-card p-6 shadow-sm transition-all hover:shadow-md">
            <div className="mb-4 flex h-10 w-10 items-center justify-center rounded-full bg-primary/10 text-primary">
              <CheckCircle className="h-5 w-5" />
            </div>
            <h3 className="text-lg font-semibold">高精度</h3>
            <p className="mt-2 text-muted-foreground">
              提供高质量的细胞类型注释，精确度可与专家手动注释相媲美。
            </p>
          </div>
        </div>
      </div>

      <div className="container mx-auto px-4 py-12">
        <h2 className="text-2xl font-bold mb-6 text-center">资源</h2>
        <div className="max-w-md mx-auto">
          <GitHubCard />
        </div>
      </div>
    </div>
  )
}
