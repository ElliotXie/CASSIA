"use client"

import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table"

const data = [
  {
    metric: "Execution Time (ms)",
    datasets: [
      { size: "10K", QuickSort: 12, MergeSort: 15, HeapSort: 18 },
      { size: "100K", QuickSort: 35, MergeSort: 42, HeapSort: 48 },
      { size: "500K", QuickSort: 75, MergeSort: 90, HeapSort: 105 },
      { size: "1M", QuickSort: 120, MergeSort: 145, HeapSort: 160 },
      { size: "5M", QuickSort: 250, MergeSort: 310, HeapSort: 340 },
      { size: "10M", QuickSort: 380, MergeSort: 470, HeapSort: 520 },
    ],
  },
  {
    metric: "Memory Usage (MB)",
    datasets: [
      { size: "10K", QuickSort: 8, MergeSort: 12, HeapSort: 10 },
      { size: "100K", QuickSort: 18, MergeSort: 25, HeapSort: 20 },
      { size: "500K", QuickSort: 30, MergeSort: 40, HeapSort: 35 },
      { size: "1M", QuickSort: 45, MergeSort: 52, HeapSort: 48 },
      { size: "5M", QuickSort: 85, MergeSort: 110, HeapSort: 95 },
      { size: "10M", QuickSort: 120, MergeSort: 160, HeapSort: 140 },
    ],
  },
  {
    metric: "Throughput (ops/sec)",
    datasets: [
      { size: "10K", QuickSort: 83000, MergeSort: 66000, HeapSort: 55000 },
      { size: "100K", QuickSort: 28000, MergeSort: 23000, HeapSort: 20000 },
      { size: "500K", QuickSort: 13000, MergeSort: 11000, HeapSort: 9500 },
      { size: "1M", QuickSort: 8300, MergeSort: 6900, HeapSort: 6250 },
      { size: "5M", QuickSort: 4000, MergeSort: 3200, HeapSort: 2900 },
      { size: "10M", QuickSort: 2600, MergeSort: 2100, HeapSort: 1900 },
    ],
  },
  {
    metric: "Time Complexity",
    value: {
      QuickSort: "O(n log n) average, O(n²) worst",
      MergeSort: "O(n log n)",
      HeapSort: "O(n log n)",
    },
  },
  {
    metric: "Space Complexity",
    value: {
      QuickSort: "O(log n)",
      MergeSort: "O(n)",
      HeapSort: "O(1)",
    },
  },
  {
    metric: "Stability",
    value: {
      QuickSort: "Not stable",
      MergeSort: "Stable",
      HeapSort: "Not stable",
    },
  },
]

export function ComparisonTable() {
  return (
    <div className="rounded-md border">
      <Table>
        <TableHeader>
          <TableRow>
            <TableHead className="w-[200px]">Metric</TableHead>
            <TableHead>Dataset Size</TableHead>
            <TableHead>QuickSort</TableHead>
            <TableHead>MergeSort</TableHead>
            <TableHead>HeapSort</TableHead>
          </TableRow>
        </TableHeader>
        <TableBody>
          {data.map((item, index) =>
            item.datasets ? (
              item.datasets.map((dataset, datasetIndex) => (
                <TableRow key={`${index}-${datasetIndex}`}>
                  {datasetIndex === 0 ? (
                    <TableCell rowSpan={item.datasets.length} className="font-medium align-top">
                      {item.metric}
                    </TableCell>
                  ) : null}
                  <TableCell>{dataset.size}</TableCell>
                  <TableCell>{dataset.QuickSort.toLocaleString()}</TableCell>
                  <TableCell>{dataset.MergeSort.toLocaleString()}</TableCell>
                  <TableCell>{dataset.HeapSort.toLocaleString()}</TableCell>
                </TableRow>
              ))
            ) : (
              <TableRow key={index}>
                <TableCell className="font-medium">{item.metric}</TableCell>
                <TableCell>-</TableCell>
                <TableCell>{item.value.QuickSort}</TableCell>
                <TableCell>{item.value.MergeSort}</TableCell>
                <TableCell>{item.value.HeapSort}</TableCell>
              </TableRow>
            ),
          )}
        </TableBody>
      </Table>
    </div>
  )
}
