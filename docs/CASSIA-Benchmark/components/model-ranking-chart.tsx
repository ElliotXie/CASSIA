"use client"

import { CartesianGrid, Legend, ResponsiveContainer, Scatter, ScatterChart, Tooltip, XAxis, YAxis, ZAxis, Label, Text } from "recharts"

// Combined data with performance and cost metrics
const modelRankingData = [
  {
    model: "Llama 4 Maverick",
    outputCost: 0.60,
    averageScore: 86.68,
    costRank: 2, // Tied with Gemini 2.5 flash
    scoreRank: 1
  },
  {
    model: "Gemini 2.5 flash",
    outputCost: 0.60,
    averageScore: 86.40,
    costRank: 2, // Tied with Llama 4 Maverick
    scoreRank: 2
  },
  {
    model: "Claude 3.7",
    outputCost: 15.00,
    averageScore: 85.97,
    costRank: 8,
    scoreRank: 3
  },
  {
    model: "Gemini 2.5 pro",
    outputCost: 10.00,
    averageScore: 85.32,
    costRank: 7,
    scoreRank: 4
  },
  {
    model: "Deepseek v3",
    outputCost: 1.10,
    averageScore: 85.27,
    costRank: 3,
    scoreRank: 5
  },
  {
    model: "GPT-4.1",
    outputCost: 8.00,
    averageScore: 84.18,
    costRank: 6,
    scoreRank: 6
  },
  {
    model: "GPT-O4 Mini High",
    outputCost: 4.40,
    averageScore: 84.14,
    costRank: 5,
    scoreRank: 7
  },
  {
    model: "QWEN3-235b",
    outputCost: 2.00,
    averageScore: 83.82,
    costRank: 4,
    scoreRank: 8
  }
]

// Color mapping for models (same as other charts)
const modelColors = {
  "Llama 4 Maverick": "#6366f1", // indigo
  "Claude 3.7": "#8b5cf6", // violet
  "Deepseek v3": "#ec4899", // pink
  "Gemini 2.5 flash": "#fbbf24", // lighter amber
  "Gemini 2.5 pro": "#f59e0b", // amber
  "GPT-4.1": "#0ea5e9", // sky blue
  "GPT-O4 Mini High": "#60a5fa", // lighter blue
  "QWEN3-235b": "#14b8a6", // teal
}

// Custom label component for each point
const CustomizedLabel = (props: any) => {
  const { cx, cy, name, fill } = props;
  
  // Get abbreviated model name for label to avoid overlap
  let shortName = name;
  if (name === "Llama 4 Maverick") shortName = "Llama 4";
  if (name === "Gemini 2.5 flash") shortName = "Gemini Flash";
  if (name === "Gemini 2.5 pro") shortName = "Gemini Pro";
  if (name === "GPT-O4 Mini High") shortName = "O4 Mini";
  if (name === "QWEN3-235b") shortName = "QWEN3";
  
  return (
    <text 
      x={cx + 10} 
      y={cy} 
      dy={-10} 
      fill={fill} 
      fontSize={11}
      fontWeight={500} 
      textAnchor="start"
    >
      {shortName}
    </text>
  );
};

// Custom tooltip component that only shows the data we want
const CustomTooltip = ({ active, payload }: any) => {
  if (active && payload && payload.length) {
    const data = payload[0].payload;
    return (
      <div className="bg-white p-3 rounded-md shadow-md border border-gray-200">
        <p className="font-semibold">{data.model}</p>
        <p className="text-sm">Average Score: <span className="font-medium">{data.averageScore.toFixed(2)}%</span></p>
        <p className="text-sm">Output Cost: <span className="font-medium">${data.outputCost.toFixed(2)}</span> per 1M tokens</p>
      </div>
    );
  }
  return null;
};

export function ModelRankingChart() {
  // Calculate the size of each point based on its overall ranking (score rank + cost rank)
  // Lower number = better overall ranking
  const sizeData = modelRankingData.map(item => ({
    ...item,
    overallRanking: item.costRank + item.scoreRank,
    size: 20 - ((item.costRank + item.scoreRank) / 2) * 2 // Size decreases as overall ranking number increases
  }))
  
  return (
    <div className="w-full h-[500px]">
      <ResponsiveContainer width="100%" height="100%">
        <ScatterChart
          margin={{ top: 20, right: 30, left: 20, bottom: 100 }} // Increased bottom margin
        >
          <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
          <XAxis 
            type="number"
            dataKey="outputCost" 
            name="Output Cost" 
            domain={[0, 16]}
            label={{ 
              value: 'Output Cost ($ per 1M tokens)', 
              position: 'insideBottom', 
              offset: -10,
              style: { textAnchor: 'middle', fill: '#4b5563', fontWeight: 500 }
            }}
            tick={{ fill: "#4b5563" }}
          />
          <YAxis 
            type="number"
            dataKey="averageScore" 
            name="Average Score" 
            domain={[83, 88]}
            label={{
              value: "Average Score (%)",
              angle: -90,
              position: "insideLeft",
              offset: 10,
              style: { fill: "#4b5563", textAnchor: "middle", fontWeight: 500 },
            }}
            tick={{ fill: "#4b5563" }}
          />
          <ZAxis 
            type="number"
            dataKey="size" 
            range={[60, 200]} 
          />
          <Tooltip content={<CustomTooltip />} cursor={{ strokeDasharray: '3 3' }} />
          <Legend 
            formatter={(value, entry) => {
              const model = sizeData.find(item => item.model === value);
              return <span style={{ color: "#4b5563", fontWeight: 500 }}>{value}</span>;
            }}
            verticalAlign="bottom"
            height={36}
            wrapperStyle={{ paddingTop: "15px" }}
          />
          {sizeData.map((item) => (
            <Scatter 
              key={item.model} 
              name={item.model} 
              data={[item]} 
              fill={modelColors[item.model as keyof typeof modelColors]}
              line={{ stroke: modelColors[item.model as keyof typeof modelColors], strokeWidth: 1 }}
              shape="circle"
              label={<CustomizedLabel fill={modelColors[item.model as keyof typeof modelColors]} />}
            />
          ))}
        </ScatterChart>
      </ResponsiveContainer>
      <div className="text-center text-sm text-gray-500 mt-4">
        <p>Optimal models appear in the top-left corner (higher score, lower cost). Larger circles indicate better overall ranking.</p>
      </div>
    </div>
  )
} 