"use client"

import { Bar, BarChart, CartesianGrid, Cell, Legend, ResponsiveContainer, Tooltip, XAxis, YAxis } from "recharts"

// Updated model data for both methods
const modelData = {
  gptcelltype: {
    "Human Kidney": [
      { model: "GPT-4o", score: 0.741 },
      { model: "GPT-4", score: 0.481 },
    ],
    "Human Lung": [
      { model: "GPT-4o", score: 0.733 },
      { model: "GPT-4", score: 0.800 },
    ],
    "Human Large Intestine": [
      { model: "GPT-4o", score: 0.294 },
      { model: "GPT-4", score: 0.324 },
    ],
    "Human Fetal Skin": [
      { model: "GPT-4o", score: 0.818 },
      { model: "GPT-4", score: 0.727 },
    ],
    "Mouse Atlas": [
      { model: "GPT-4o", score: 0.720 },
      { model: "GPT-4", score: 0.700 },
    ],
  },
  cassia: {
    "Human Kidney": [
      { model: "Llama 4 Maverick", score: 0.8704 },
      { model: "Gemini 2.5 flash", score: 0.8389 },
      { model: "Claude 3.7", score: 0.8667 },
      { model: "Gemini 2.5 pro", score: 0.8222 },
      { model: "Deepseek v3", score: 0.8444 },
      { model: "GPT-4.1", score: 0.8167 },
      { model: "GPT-O4 Mini High", score: 0.8148 },
      { model: "QWEN3-235b", score: 0.8278 },
    ],
    "Human Lung": [
      { model: "Llama 4 Maverick", score: 0.94 },
      { model: "Gemini 2.5 flash", score: 0.9267 },
      { model: "Deepseek v3", score: 0.9133 },
      { model: "Claude 3.7", score: 0.9033 },
      { model: "Gemini 2.5 pro", score: 0.9 },
      { model: "GPT-4.1", score: 0.9 },
      { model: "GPT-O4 Mini High", score: 0.9 },
      { model: "QWEN3-235b", score: 0.8833 },
    ],
    "Human Large Intestine": [
      { model: "Gemini 2.5 flash", score: 0.9 },
      { model: "Gemini 2.5 pro", score: 0.8938 },
      { model: "GPT-4.1", score: 0.875 },
      { model: "Deepseek v3", score: 0.8688 },
      { model: "QWEN3-235b", score: 0.8594 },
      { model: "Llama 4 Maverick", score: 0.85 },
      { model: "Claude 3.7", score: 0.8312 },
      { model: "GPT-O4 Mini High", score: 0.825 },
    ],
    "Human Fetal Skin": [
      { model: "GPT-O4 Mini High", score: 0.8591 },
      { model: "Claude 3.7", score: 0.8545 },
      { model: "Llama 4 Maverick", score: 0.8455 },
      { model: "Gemini 2.5 flash", score: 0.8273 },
      { model: "QWEN3-235b", score: 0.8273 },
      { model: "Deepseek v3", score: 0.8182 },
      { model: "GPT-4.1", score: 0.8182 },
      { model: "Gemini 2.5 pro", score: 0.8182 },
    ],
    "Mouse Atlas": [
      { model: "Gemini 2.5 pro", score: 0.8479 },
      { model: "Gemini 2.5 flash", score: 0.8458 },
      { model: "Claude 3.7", score: 0.8458 },
      { model: "GPT-O4 Mini High", score: 0.8375 },
      { model: "Llama 4 Maverick", score: 0.8375 },
      { model: "Deepseek v3", score: 0.8292 },
      { model: "GPT-4.1", score: 0.825 },
      { model: "QWEN3-235b", score: 0.8125 },
    ],
  },
}

// Updated color mapping for models
const modelColors = {
  // GPTcelltype models
  "GPT-4o": "#3b82f6", // blue
  "GPT-4": "#10b981", // green

  // CASSIA models
  "GPT-4.1": "#0ea5e9", // sky blue
  "GPT-O4 Mini High": "#60a5fa", // lighter blue
  "Claude 3.7": "#8b5cf6", // violet
  "Gemini 2.5 pro": "#f59e0b", // amber
  "Gemini 2.5 flash": "#fbbf24", // lighter amber
  "Deepseek v3": "#ec4899", // pink
  "Llama 4 Maverick": "#6366f1", // indigo
  "QWEN3-235b": "#14b8a6", // teal
}

interface ModelComparisonChartProps {
  method: "gptcelltype" | "cassia"
  selectedTissue: string
}

export function ModelComparisonChart({ method, selectedTissue }: ModelComparisonChartProps) {
  // Get data for the selected tissue or combine all tissues
  const getChartData = () => {
    if (selectedTissue === "All Tissues") {
      // For "All Tissues", we'll show data for all tissues
      const tissues = ["Human Kidney", "Human Lung", "Human Large Intestine", "Human Fetal Skin", "Mouse Atlas"]

      return tissues.map((tissue) => {
        // Find the data points for this tissue
        const tissueData = modelData[method][tissue]
        const dataPoint: any = { tissue }

        // Add each model's score for this tissue
        tissueData.forEach((item) => {
          dataPoint[item.model] = item.score
        })

        return dataPoint
      })
    } else {
      // For a specific tissue, transform the data for the chart
      return modelData[method][selectedTissue].map((item) => ({
        model: item.model,
        score: item.score,
        color: modelColors[item.model as keyof typeof modelColors],
      }))
      // Sort by score in descending order
      .sort((a: { score: number }, b: { score: number }) => b.score - a.score)
    }
  }

  // Get the chart data based on the selected tissue
  const chartData = getChartData()

  // Get the list of models for the current method
  const getModelsForMethod = () => {
    // Return models in the order of the average score rank
    if (method === "cassia") {
      // This order matches the scoreRank in the Cost Performance Analysis chart
      // Ensures consistent ordering between charts while preserving tissue-specific ordering
      return [
        "Llama 4 Maverick",
        "Gemini 2.5 flash",
        "Claude 3.7",
        "Gemini 2.5 pro", 
        "Deepseek v3",
        "GPT-4.1",
        "GPT-O4 Mini High",
        "QWEN3-235b"
      ];
    }
    // For gptcelltype, use the original ordering
    return modelData[method]["Human Kidney"].map((item) => item.model)
  }

  const methodModels = getModelsForMethod()

  return (
    <div className="w-full h-[400px]">
      <ResponsiveContainer width="100%" height="100%">
        {selectedTissue === "All Tissues" ? (
          // Show bar chart with all tissues on x-axis when "All Tissues" is selected
          <BarChart
            data={chartData}
            margin={{
              top: 20,
              right: 30,
              left: 20,
              bottom: 70, // Increased bottom margin to accommodate horizontal labels
            }}
          >
            <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
            <XAxis
              dataKey="tissue"
              angle={0} // Changed from -45 to 0 for horizontal labels
              textAnchor="middle"
              height={80}
              tick={{ fill: "#4b5563" }}
              axisLine={{ stroke: "#e5e7eb" }}
              tickMargin={10} // Added margin for better spacing
            />
            <YAxis
              domain={[0, 1]}
              label={{
                value: "Accuracy Score",
                angle: -90,
                position: "insideLeft",
                style: { fill: "#4b5563", textAnchor: "middle" },
              }}
              tick={{ fill: "#4b5563" }}
              axisLine={{ stroke: "#e5e7eb" }}
              tickLine={{ stroke: "#e5e7eb" }}
            />
            <Tooltip
              formatter={(value, name) => [typeof value === 'number' ? value.toFixed(2) : value, name]}
              contentStyle={{
                backgroundColor: "white",
                borderRadius: "0.5rem",
                boxShadow: "0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)",
                border: "1px solid #e5e7eb",
              }}
              labelStyle={{ fontWeight: "bold", marginBottom: "0.5rem" }}
              itemStyle={{ padding: "2px 0" }}
              isAnimationActive={false}
            />
            <Legend
              wrapperStyle={{ paddingTop: "20px" }}
              formatter={(value) => <span style={{ color: "#4b5563", fontWeight: 500 }}>{value}</span>}
            />
            {methodModels.map((model) => (
              <Bar
                key={model}
                dataKey={model}
                name={model}
                fill={modelColors[model as keyof typeof modelColors]}
                radius={[4, 4, 0, 0]}
                animationDuration={1500}
              />
            ))}
          </BarChart>
        ) : (
          // Show horizontal bar chart for a specific tissue
          <BarChart
            data={chartData}
            margin={{
              top: 20,
              right: 30,
              left: 150, // Increased left margin for longer model names
              bottom: 20,
            }}
            layout="vertical"
          >
            <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" horizontal={false} />
            <XAxis
              type="number"
              domain={[0, 1]}
              tick={{ fill: "#4b5563" }}
              axisLine={{ stroke: "#e5e7eb" }}
              tickLine={{ stroke: "#e5e7eb" }}
            />
            <YAxis
              dataKey="model"
              type="category"
              width={150} // Increased width for longer model names
              tick={{ fill: "#4b5563" }}
              axisLine={{ stroke: "#e5e7eb" }}
            />
            <Tooltip
              formatter={(value) => [typeof value === 'number' ? value.toFixed(2) : value, "Accuracy Score"]}
              contentStyle={{
                backgroundColor: "white",
                borderRadius: "0.5rem",
                boxShadow: "0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)",
                border: "1px solid #e5e7eb",
              }}
              labelStyle={{ fontWeight: "bold", marginBottom: "0.5rem" }}
              itemStyle={{ padding: "2px 0" }}
              isAnimationActive={false}
            />
            <Legend
              wrapperStyle={{ paddingTop: "20px" }}
              formatter={(value) => <span style={{ color: "#4b5563", fontWeight: 500 }}>{value}</span>}
            />
            <Bar
              dataKey="score"
              name="Accuracy Score"
              radius={[0, 4, 4, 0]}
              animationDuration={1500}
              fill="#3b82f6"
              // Use the model-specific colors
              isAnimationActive={true}
              fillOpacity={0.9}
            >
              {chartData.map((entry: { model: string; score: number; color: string }, index: number) => (
                <Cell key={`cell-${index}`} fill={entry.color} />
              ))}
            </Bar>
          </BarChart>
        )}
      </ResponsiveContainer>
    </div>
  )
}
