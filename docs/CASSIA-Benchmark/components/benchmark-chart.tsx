"use client"

import { Bar, BarChart, CartesianGrid, Legend, ResponsiveContainer, Tooltip, XAxis, YAxis } from "recharts"

const data = [
  { name: "QuickSort", executionTime: 120, memoryUsage: 45, throughput: 8300 },
  { name: "MergeSort", executionTime: 145, memoryUsage: 52, throughput: 6900 },
  { name: "HeapSort", executionTime: 160, memoryUsage: 48, throughput: 6250 },
  { name: "BubbleSort", executionTime: 450, memoryUsage: 40, throughput: 2200 },
  { name: "InsertionSort", executionTime: 380, memoryUsage: 38, throughput: 2600 },
  { name: "RadixSort", executionTime: 110, memoryUsage: 65, throughput: 9100 },
  { name: "CountingSort", executionTime: 95, memoryUsage: 72, throughput: 10500 },
  { name: "BucketSort", executionTime: 135, memoryUsage: 58, throughput: 7400 },
]

export function BenchmarkChart() {
  return (
    <div className="w-full h-[400px]">
      <ResponsiveContainer width="100%" height="100%">
        <BarChart
          data={data}
          margin={{
            top: 20,
            right: 30,
            left: 20,
            bottom: 60,
          }}
        >
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis dataKey="name" angle={-45} textAnchor="end" height={60} />
          <YAxis />
          <Tooltip />
          <Legend />
          <Bar dataKey="executionTime" name="Execution Time (ms)" fill="#3b82f6" />
        </BarChart>
      </ResponsiveContainer>
    </div>
  )
}
