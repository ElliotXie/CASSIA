import Link from "next/link"
import { Button } from "@/components/ui/button"
import { Dna } from "lucide-react"

export default function NotFound() {
  return (
    <div className="flex min-h-screen flex-col items-center justify-center p-4 text-center">
      <div className="flex h-16 w-16 items-center justify-center rounded-full bg-primary/10 text-primary mb-6">
        <Dna className="h-8 w-8" />
      </div>
      <h1 className="text-4xl font-bold tracking-tight sm:text-5xl">404 - Page Not Found</h1>
      <p className="mt-4 text-lg text-muted-foreground">The page you are looking for does not exist.</p>
      <Button asChild className="mt-8">
        <Link href="/">Return Home</Link>
      </Button>
    </div>
  )
}
