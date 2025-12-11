import Link from "next/link"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardDescription, CardFooter, CardHeader, CardTitle } from "@/components/ui/card"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Textarea } from "@/components/ui/textarea"
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"
import { ArrowLeft } from "lucide-react"

export default function SubmitBenchmark() {
  return (
    <div className="container mx-auto py-10 max-w-2xl">
      <div className="mb-8">
        <Link
          href="/"
          className="flex items-center text-muted-foreground hover:text-foreground transition-colors gap-1 mb-4"
        >
          <ArrowLeft className="h-4 w-4" />
          Back to Dashboard
        </Link>
        <h1 className="text-3xl font-bold tracking-tight">Submit Benchmark Results</h1>
        <p className="text-muted-foreground mt-2">Add your benchmark data to compare with other methods</p>
      </div>

      <Card>
        <CardHeader>
          <CardTitle>Benchmark Details</CardTitle>
          <CardDescription>Please provide information about your benchmark method and results</CardDescription>
        </CardHeader>
        <CardContent className="space-y-6">
          <div className="space-y-2">
            <Label htmlFor="method-name">Method Name</Label>
            <Input id="method-name" placeholder="e.g., Quick Sort, Binary Search, etc." />
          </div>

          <div className="space-y-2">
            <Label htmlFor="category">Category</Label>
            <Select>
              <SelectTrigger>
                <SelectValue placeholder="Select category" />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="sorting">Sorting Algorithms</SelectItem>
                <SelectItem value="search">Search Algorithms</SelectItem>
                <SelectItem value="graph">Graph Algorithms</SelectItem>
                <SelectItem value="ml">Machine Learning</SelectItem>
                <SelectItem value="crypto">Cryptography</SelectItem>
                <SelectItem value="other">Other</SelectItem>
              </SelectContent>
            </Select>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label htmlFor="execution-time">Execution Time (ms)</Label>
              <Input id="execution-time" type="number" placeholder="e.g., 125" />
            </div>
            <div className="space-y-2">
              <Label htmlFor="memory-usage">Memory Usage (MB)</Label>
              <Input id="memory-usage" type="number" placeholder="e.g., 42" />
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label htmlFor="throughput">Throughput (ops/sec)</Label>
              <Input id="throughput" type="number" placeholder="e.g., 10000" />
            </div>
            <div className="space-y-2">
              <Label htmlFor="dataset-size">Dataset Size</Label>
              <Input id="dataset-size" type="number" placeholder="e.g., 1000000" />
            </div>
          </div>

          <div className="space-y-2">
            <Label htmlFor="environment">Environment</Label>
            <Input id="environment" placeholder="e.g., Node.js v18.0.0, Intel i7-12700K, 32GB RAM" />
          </div>

          <div className="space-y-2">
            <Label htmlFor="description">Description & Implementation Details</Label>
            <Textarea
              id="description"
              placeholder="Describe your implementation, optimizations, and any other relevant details"
              rows={5}
            />
          </div>

          <div className="space-y-2">
            <Label htmlFor="code-snippet">Code Snippet (optional)</Label>
            <Textarea id="code-snippet" placeholder="Paste your code here" rows={8} className="font-mono text-sm" />
          </div>
        </CardContent>
        <CardFooter className="flex justify-end gap-2">
          <Button variant="outline">Cancel</Button>
          <Button>Submit Benchmark</Button>
        </CardFooter>
      </Card>
    </div>
  )
}
