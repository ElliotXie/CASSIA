import Link from "next/link"
import { FileBarChart } from "lucide-react"

export function Navbar() {
  return (
    <header className="navbar sticky top-0 z-10 border-b">
      <div className="container mx-auto flex h-16 items-center justify-between">
        <Link href="/" className="flex items-center gap-2">
          <FileBarChart className="h-6 w-6 text-blue-600" />
          <span className="font-bold text-xl bg-clip-text text-transparent bg-gradient-to-r from-blue-600 to-emerald-600">
            LLM Cell Annotation Benchmark
          </span>
        </Link>
        <nav className="flex items-center gap-8">
          <Link href="/" className="text-sm font-medium text-gray-700 hover:text-blue-600 transition-colors">
            Home
          </Link>
          <Link
            href="/methods/cassia"
            className="text-sm font-medium text-gray-700 hover:text-emerald-600 transition-colors"
          >
            CASSIA
          </Link>
          <Link
            href="/methods/gptcelltype"
            className="text-sm font-medium text-gray-700 hover:text-blue-600 transition-colors"
          >
            GPTcelltype
          </Link>
        </nav>
      </div>
    </header>
  )
}
