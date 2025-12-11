"use client"

import {
  Bar,
  BarChart,
  CartesianGrid,
  Legend,
  Line,
  LineChart,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts"

const executionTimeData = [
  { name: "10K", QuickSort: 12, MergeSort: 15, HeapSort: 18 },
  { name: "100K", QuickSort: 35, MergeSort: 42, HeapSort: 48 },
  { name: "500K", QuickSort: 75, MergeSort: 90, HeapSort: 105 },
  { name: "1M", QuickSort: 120, MergeSort: 145, HeapSort: 160 },
  { name: "5M", QuickSort: 250, MergeSort: 310, HeapSort: 340 },
  { name: "10M", QuickSort: 380, MergeSort: 470, HeapSort: 520 },
]

const memoryUsageData = [
  { name: "10K", QuickSort: 8, MergeSort: 12, HeapSort: 10 },
  { name: "100K", QuickSort: 18, MergeSort: 25, HeapSort: 20 },
  { name: "500K", QuickSort: 30, MergeSort: 40, HeapSort: 35 },
  { name: "1M", QuickSort: 45, MergeSort: 52, HeapSort: 48 },
  { name: "5M", QuickSort: 85, MergeSort: 110, HeapSort: 95 },
  { name: "10M", QuickSort: 120, MergeSort: 160, HeapSort: 140 },
]

const throughputData = [
  { name: "10K", QuickSort: 83000, MergeSort: 66000, HeapSort: 55000 },
  { name: "100K", QuickSort: 28000, MergeSort: 23000, HeapSort: 20000 },
  { name: "500K", QuickSort: 13000, MergeSort: 11000, HeapSort: 9500 },
  { name: "1M", QuickSort: 8300, MergeSort: 6900, HeapSort: 6250 },
  { name: "5M", QuickSort: 4000, MergeSort: 3200, HeapSort: 2900 },
  { name: "10M", QuickSort: 2600, MergeSort: 2100, HeapSort: 1900 },
]

const scalingData = [
  { name: "10K", QuickSort: 1, MergeSort: 1, HeapSort: 1 },
  { name: "100K", QuickSort: 2.9, MergeSort: 2.8, HeapSort: 2.7 },
  { name: "500K", QuickSort: 6.3, MergeSort: 6, HeapSort: 5.8 },
  { name: "1M", QuickSort: 10, MergeSort: 9.7, HeapSort: 8.9 },
  { name: "5M", QuickSort: 20.8, MergeSort: 20.7, HeapSort: 18.9 },
  { name: "10M", QuickSort: 31.7, MergeSort: 31.3, HeapSort: 28.9 },
]

type MetricType = "execution-time" | "memory-usage" | "throughput" | "scaling"

interface ComparisonChartProps {
  metric: MetricType
}

export function ComparisonChart({ metric }: ComparisonChartProps) {
  const getChartData = () => {
    switch (metric) {
      case "execution-time":
        return {
          data: executionTimeData,
          yAxisLabel: "Time (ms)",
          title: "Execution Time Comparison",
        }
      case "memory-usage":
        return {
          data: memoryUsageData,
          yAxisLabel: "Memory (MB)",
          title: "Memory Usage Comparison",
        }
      case "throughput":
        return {
          data: throughputData,
          yAxisLabel: "Operations/sec",
          title: "Throughput Comparison",
        }
      case "scaling":
        return {
          data: scalingData,
          yAxisLabel: "Relative Performance",
          title: "Scaling Comparison",
        }
      default:
        return {
          data: executionTimeData,
          yAxisLabel: "Time (ms)",
          title: "Execution Time Comparison",
        }
    }
  }

  const { data, yAxisLabel, title } = getChartData()

  // Use LineChart for scaling, BarChart for others
  return (
    <div className="w-full">
      <h3 className="text-lg font-medium mb-4">{title}</h3>
      <div className="h-[400px]">
        <ResponsiveContainer width="100%" height="100%">
          {metric === "scaling" ? (
            <LineChart
              data={data}
              margin={{
                top: 20,
                right: 30,
                left: 20,
                bottom: 20,
              }}
            >
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="name" />
              <YAxis label={{ value: yAxisLabel, angle: -90, position: "insideLeft" }} />
              <Tooltip />
              <Legend />
              <Line type="monotone" dataKey="QuickSort" stroke="#3b82f6" strokeWidth={2} />
              <Line type="monotone" dataKey="MergeSort" stroke="#10b981" strokeWidth={2} />
              <Line type="monotone" dataKey="HeapSort" stroke="#f59e0b" strokeWidth={2} />
            </LineChart>
          ) : (
            <BarChart
              data={data}
              margin={{
                top: 20,
                right: 30,
                left: 20,
                bottom: 20,
              }}
            >
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="name" />
              <YAxis label={{ value: yAxisLabel, angle: -90, position: "insideLeft" }} />
              <Tooltip />
              <Legend />
              <Bar dataKey="QuickSort" fill="#3b82f6" />
              <Bar dataKey="MergeSort" fill="#10b981" />
              <Bar dataKey="HeapSort" fill="#f59e0b" />
            </BarChart>
          )}
        </ResponsiveContainer>
      </div>
    </div>
  )
}
