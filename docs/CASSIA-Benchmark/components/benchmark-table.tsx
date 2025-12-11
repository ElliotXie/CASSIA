"use client"

import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table"
import { Badge } from "@/components/ui/badge"
import { ArrowUpDown, ArrowUp, ArrowDown } from "lucide-react"
import { Button } from "@/components/ui/button"
import { useState } from "react"

const data = [
  {
    id: 1,
    name: "QuickSort",
    category: "Sorting",
    executionTime: 120,
    memoryUsage: 45,
    throughput: 8300,
    datasetSize: "1M",
    environment: "Node.js v18",
    date: "2023-04-15",
  },
  {
    id: 2,
    name: "MergeSort",
    category: "Sorting",
    executionTime: 145,
    memoryUsage: 52,
    throughput: 6900,
    datasetSize: "1M",
    environment: "Node.js v18",
    date: "2023-04-14",
  },
  {
    id: 3,
    name: "HeapSort",
    category: "Sorting",
    executionTime: 160,
    memoryUsage: 48,
    throughput: 6250,
    datasetSize: "1M",
    environment: "Node.js v18",
    date: "2023-04-12",
  },
  {
    id: 4,
    name: "BubbleSort",
    category: "Sorting",
    executionTime: 450,
    memoryUsage: 40,
    throughput: 2200,
    datasetSize: "1M",
    environment: "Node.js v18",
    date: "2023-04-10",
  },
  {
    id: 5,
    name: "RadixSort",
    category: "Sorting",
    executionTime: 110,
    memoryUsage: 65,
    throughput: 9100,
    datasetSize: "1M",
    environment: "Node.js v18",
    date: "2023-04-18",
  },
  {
    id: 6,
    name: "CountingSort",
    category: "Sorting",
    executionTime: 95,
    memoryUsage: 72,
    throughput: 10500,
    datasetSize: "1M",
    environment: "Node.js v18",
    date: "2023-04-20",
  },
]

type SortDirection = "asc" | "desc" | null
type SortField = "name" | "executionTime" | "memoryUsage" | "throughput" | "date" | null

export function BenchmarkTable() {
  const [sortField, setSortField] = useState<SortField>("executionTime")
  const [sortDirection, setSortDirection] = useState<SortDirection>("asc")

  const handleSort = (field: SortField) => {
    if (sortField === field) {
      setSortDirection(sortDirection === "asc" ? "desc" : sortDirection === "desc" ? null : "asc")
      if (sortDirection === null) {
        setSortField(null)
      }
    } else {
      setSortField(field)
      setSortDirection("asc")
    }
  }

  const sortedData = [...data].sort((a, b) => {
    if (sortField === null || sortDirection === null) return 0

    const aValue = a[sortField]
    const bValue = b[sortField]

    if (sortDirection === "asc") {
      return aValue > bValue ? 1 : -1
    } else {
      return aValue < bValue ? 1 : -1
    }
  })

  const getSortIcon = (field: SortField) => {
    if (sortField !== field) return <ArrowUpDown className="h-4 w-4" />
    if (sortDirection === "asc") return <ArrowUp className="h-4 w-4" />
    if (sortDirection === "desc") return <ArrowDown className="h-4 w-4" />
    return <ArrowUpDown className="h-4 w-4" />
  }

  return (
    <div className="rounded-md border">
      <Table>
        <TableHeader>
          <TableRow>
            <TableHead>
              <Button variant="ghost" size="sm" onClick={() => handleSort("name")} className="gap-1 font-medium">
                Method {getSortIcon("name")}
              </Button>
            </TableHead>
            <TableHead>Category</TableHead>
            <TableHead>
              <Button
                variant="ghost"
                size="sm"
                onClick={() => handleSort("executionTime")}
                className="gap-1 font-medium"
              >
                Execution Time (ms) {getSortIcon("executionTime")}
              </Button>
            </TableHead>
            <TableHead>
              <Button variant="ghost" size="sm" onClick={() => handleSort("memoryUsage")} className="gap-1 font-medium">
                Memory (MB) {getSortIcon("memoryUsage")}
              </Button>
            </TableHead>
            <TableHead>
              <Button variant="ghost" size="sm" onClick={() => handleSort("throughput")} className="gap-1 font-medium">
                Throughput (ops/sec) {getSortIcon("throughput")}
              </Button>
            </TableHead>
            <TableHead>Dataset</TableHead>
            <TableHead>Environment</TableHead>
            <TableHead>
              <Button variant="ghost" size="sm" onClick={() => handleSort("date")} className="gap-1 font-medium">
                Date {getSortIcon("date")}
              </Button>
            </TableHead>
          </TableRow>
        </TableHeader>
        <TableBody>
          {sortedData.map((item) => (
            <TableRow key={item.id}>
              <TableCell className="font-medium">{item.name}</TableCell>
              <TableCell>
                <Badge variant="outline">{item.category}</Badge>
              </TableCell>
              <TableCell>{item.executionTime}</TableCell>
              <TableCell>{item.memoryUsage}</TableCell>
              <TableCell>{item.throughput}</TableCell>
              <TableCell>{item.datasetSize}</TableCell>
              <TableCell>{item.environment}</TableCell>
              <TableCell>{item.date}</TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  )
}
