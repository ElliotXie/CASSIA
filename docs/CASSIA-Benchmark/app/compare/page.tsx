import Link from "next/link"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import { ComparisonChart } from "@/components/comparison-chart"
import { ComparisonTable } from "@/components/comparison-table"
import { ArrowLeft, Download } from "lucide-react"

export default function CompareBenchmarks() {
  return (
    <div className="container mx-auto py-10 space-y-8">
      <div className="flex justify-between items-center">
        <div>
          <Link
            href="/"
            className="flex items-center text-muted-foreground hover:text-foreground transition-colors gap-1 mb-4"
          >
            <ArrowLeft className="h-4 w-4" />
            Back to Dashboard
          </Link>
          <h1 className="text-3xl font-bold tracking-tight">Compare Benchmarks</h1>
          <p className="text-muted-foreground mt-2">Side-by-side comparison of selected methods</p>
        </div>
        <Button variant="outline" className="gap-2">
          <Download className="h-4 w-4" />
          Export Data
        </Button>
      </div>

      <Card>
        <CardHeader>
          <CardTitle>Performance Comparison</CardTitle>
          <CardDescription>Comparing QuickSort, MergeSort, and HeapSort algorithms</CardDescription>
        </CardHeader>
        <CardContent>
          <Tabs defaultValue="execution-time">
            <TabsList className="mb-4">
              <TabsTrigger value="execution-time">Execution Time</TabsTrigger>
              <TabsTrigger value="memory-usage">Memory Usage</TabsTrigger>
              <TabsTrigger value="throughput">Throughput</TabsTrigger>
              <TabsTrigger value="scaling">Scaling</TabsTrigger>
            </TabsList>
            <TabsContent value="execution-time">
              <ComparisonChart metric="execution-time" />
            </TabsContent>
            <TabsContent value="memory-usage">
              <ComparisonChart metric="memory-usage" />
            </TabsContent>
            <TabsContent value="throughput">
              <ComparisonChart metric="throughput" />
            </TabsContent>
            <TabsContent value="scaling">
              <ComparisonChart metric="scaling" />
            </TabsContent>
          </Tabs>
        </CardContent>
      </Card>

      <Card>
        <CardHeader>
          <CardTitle>Detailed Metrics</CardTitle>
          <CardDescription>Complete benchmark data for all metrics</CardDescription>
        </CardHeader>
        <CardContent>
          <ComparisonTable />
        </CardContent>
      </Card>
    </div>
  )
}
