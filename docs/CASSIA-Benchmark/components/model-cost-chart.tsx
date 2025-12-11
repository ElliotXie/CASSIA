"use client"

import { Bar, BarChart, CartesianGrid, Legend, ResponsiveContainer, Tooltip, XAxis, YAxis } from "recharts"

// Cost data for models
const costData = [
  {
    model: "Llama 4 Maverick",
    inputCost: 0.17,
    outputCost: 0.60,
    imageCost: 0.67
  },
  {
    model: "Claude 3.7",
    inputCost: 3.00,
    outputCost: 15.00,
    imageCost: 4.80
  },
  {
    model: "Deepseek v3",
    inputCost: 0.27,
    outputCost: 1.10,
    imageCost: 0
  },
  {
    model: "Gemini 2.5 flash",
    inputCost: 0.15,
    outputCost: 0.60,
    imageCost: 0.62
  },
  {
    model: "Gemini 2.5 pro",
    inputCost: 1.25,
    outputCost: 10.00,
    imageCost: 5.16
  },
  {
    model: "GPT-4.1",
    inputCost: 2.00,
    outputCost: 8.00,
    imageCost: 0
  },
  {
    model: "GPT-O4 Mini High",
    inputCost: 1.10,
    outputCost: 4.40,
    imageCost: 0.84
  }
]

// Transform data for chart
const transformedCostData = [
  {
    name: "Input Cost",
    description: "per 1M tokens ($)",
    ...costData.reduce((acc, item) => {
      acc[item.model] = item.inputCost
      return acc
    }, {})
  },
  {
    name: "Output Cost",
    description: "per 1M tokens ($)",
    ...costData.reduce((acc, item) => {
      acc[item.model] = item.outputCost
      return acc
    }, {})
  },
  {
    name: "Image Cost",
    description: "per 1K input images ($)",
    ...costData.reduce((acc, item) => {
      acc[item.model] = item.imageCost
      return acc
    }, {})
  }
]

// Color mapping for models
const modelColors = {
  "Llama 4 Maverick": "#6366f1", // indigo
  "Claude 3.7": "#8b5cf6", // violet
  "Deepseek v3": "#ec4899", // pink
  "Gemini 2.5 flash": "#fbbf24", // lighter amber
  "Gemini 2.5 pro": "#f59e0b", // amber
  "GPT-4.1": "#0ea5e9", // sky blue
  "GPT-O4 Mini High": "#60a5fa", // lighter blue
}

export function ModelCostChart() {
  const modelList = costData.map(item => item.model)
  
  return (
    <div className="w-full h-[500px]">
      <ResponsiveContainer width="100%" height="100%">
        <BarChart
          data={transformedCostData}
          margin={{ top: 20, right: 30, left: 20, bottom: 70 }}
        >
          <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
          <XAxis 
            dataKey="name"
            angle={0}
            textAnchor="middle"
            height={80}
            tick={{ fill: "#4b5563" }}
            axisLine={{ stroke: "#e5e7eb" }}
            tickMargin={10}
          />
          <YAxis
            label={{
              value: "Cost ($)",
              angle: -90,
              position: "insideLeft",
              style: { fill: "#4b5563", textAnchor: "middle" },
            }}
            tick={{ fill: "#4b5563" }}
            axisLine={{ stroke: "#e5e7eb" }}
            tickLine={{ stroke: "#e5e7eb" }}
          />
          <Tooltip
            formatter={(value, name) => [`$${value.toFixed(2)}`, name]}
            contentStyle={{
              backgroundColor: "white",
              borderRadius: "0.5rem",
              boxShadow: "0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)",
              border: "1px solid #e5e7eb",
            }}
            labelStyle={{ fontWeight: "bold", marginBottom: "0.5rem" }}
            itemStyle={{ padding: "2px 0" }}
            labelFormatter={(value) => `${value} ${transformedCostData.find(item => item.name === value)?.description || ''}`}
            isAnimationActive={false}
          />
          <Legend
            wrapperStyle={{ paddingTop: "20px" }}
            formatter={(value) => <span style={{ color: "#4b5563", fontWeight: 500 }}>{value}</span>}
          />
          {modelList.map((model) => (
            <Bar
              key={model}
              dataKey={model}
              name={model}
              fill={modelColors[model]}
              radius={[4, 4, 0, 0]}
              animationDuration={1500}
            />
          ))}
        </BarChart>
      </ResponsiveContainer>
    </div>
  )
} 