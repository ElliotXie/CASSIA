#!/usr/bin/env python3
"""
Test script for the new symphonyCompare function - the upgraded cell type comparison.

This demonstrates how to use the unified function that automatically:
- Runs multiple models in parallel
- Triggers discussion rounds when there's no consensus
- Generates beautiful HTML reports
- Provides structured output with confidence scores
"""

import os
import sys

# Add CASSIA to path if needed
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Load environment variables from .env file
try:
    from dotenv import load_dotenv
    # Look for .env file in CASSIA subdirectory
    env_path = os.path.join(os.path.dirname(__file__), 'CASSIA', '.env')
    if os.path.exists(env_path):
        load_dotenv(env_path)
        print(f"✅ Loaded environment variables from {env_path}")
    else:
        # Try loading from current directory
        load_dotenv()
        print("✅ Attempted to load .env from current directory")
except ImportError:
    print("⚠️  python-dotenv not installed. Install with: pip install python-dotenv")

try:
    from CASSIA import symphonyCompare
    print("✅ Successfully imported symphonyCompare from CASSIA")
except ImportError as e:
    print(f"❌ Failed to import: {e}")
    sys.exit(1)


def demo_symphony_compare():
    """Demonstrate the symphonyCompare function with various scenarios."""
    
    print("\n" + "="*70)
    print("🎼 CASSIA Symphony Compare - Demonstration")
    print("="*70)
    
    # Check for API key
    if not os.environ.get('OPENROUTER_API_KEY'):
        print("\n⚠️  Warning: OPENROUTER_API_KEY not set")
        print("To run this demo with real API calls, set your API key:")
        print("  export OPENROUTER_API_KEY='your-key-here'")
        print("\nFor now, showing example usage...\n")
        
        # Show example code
        print("📝 Example Usage:")
        print("-" * 50)
        print("""
import CASSIA

# Basic usage - let Symphony Compare handle everything
results = CASSIA.symphonyCompare(
    tissue="peripheral blood",
    celltypes=["T cell", "B cell", "NK cell", "Monocyte"],
    marker_set="CD3, CD4, CD8, CD19, CD20, CD16, CD56, CD14",
    species="human"
)

# Access the results
print(f"Consensus: {results['consensus']}")
print(f"Confidence: {results['confidence']:.1%}")
print(f"Report saved to: {results['html_file']}")

# Advanced usage with custom settings
results = CASSIA.symphonyCompare(
    tissue="brain",
    celltypes=["Neuron", "Astrocyte", "Microglia", "Oligodendrocyte"],
    marker_set="RBFOX3, GFAP, IBA1, OLIG2, MAP2, S100B, CD11B, MBP",
    species="mouse",
    model_preset="quartet",  # Use 4 models instead of 3
    enable_discussion=True,  # Enable automatic discussion rounds
    max_discussion_rounds=3,  # Allow up to 3 discussion rounds
    consensus_threshold=0.75,  # 75% of models must agree
    output_dir="./symphony_results",
    verbose=True
)

# Work with the structured data
df = results['dataframe']
summary = results['summary']

print(f"Total rounds: {summary['total_rounds']}")
print(f"Models used: {summary['models_used']}")

# Cell type scores
for celltype, scores in summary['celltype_scores'].items():
    print(f"{celltype}: {scores['mean']:.1f} ± {scores['std']:.1f}")
""")
        print("-" * 50)
        return
    
    # If API key is available, run actual tests
    print("\n🔑 API key detected - running live demonstration")
    
    # Test Case 1: Blood cell types (should reach consensus easily)
    print("\n" + "="*60)
    print("📋 Test Case 1: Blood Cell Types")
    print("="*60)
    
    try:
        results = symphonyCompare(
            tissue="peripheral blood",
            celltypes=["T cell", "B cell", "Monocyte"],
            marker_set="CD3, CD4, CD8, CD19, CD20, CD14, CD16",
            species="human",
            model_preset="budget",  # Use budget models for testing
            enable_discussion=True,
            max_discussion_rounds=2,
            output_dir="./symphony_test_results",
            output_basename="blood_cells_test",
            verbose=True
        )
        
        print(f"\n📊 Results Summary:")
        print(f"  • Consensus: {results['consensus'] or 'No consensus reached'}")
        print(f"  • Confidence: {results['confidence']:.1%}")
        print(f"  • Total rounds: {results['summary']['total_rounds']}")
        print(f"  • CSV file: {results['csv_file']}")
        print(f"  • HTML report: {results['html_file']}")
        
    except Exception as e:
        print(f"❌ Test Case 1 failed: {e}")
        import traceback
        traceback.print_exc()
    
    # Test Case 2: Brain cell types (might need discussion)
    print("\n" + "="*60)
    print("📋 Test Case 2: Brain Cell Types")
    print("="*60)
    
    try:
        results = symphonyCompare(
            tissue="brain cortex",
            celltypes=["Neuron", "Astrocyte", "Microglia"],
            marker_set="RBFOX3, GFAP, IBA1, S100B, MAP2, P2RY12, SYP, AQP4, TMEM119, OLIG2, ALDH1L1, CD45, NEUN, CX3CR1, TUBB3, GAD1",
            species="human",
            model_preset="budget",
            enable_discussion=True,
            max_discussion_rounds=2,
            consensus_threshold=0.7,  # Lower threshold
            output_dir="./symphony_test_results",
            output_basename="brain_cells_test",
            verbose=True
        )
        
        print(f"\n📊 Results Summary:")
        print(f"  • Consensus: {results['consensus'] or 'No consensus reached'}")
        print(f"  • Confidence: {results['confidence']:.1%}")
        
        # Show score details
        print(f"\n📈 Cell Type Scores:")
        for ct, scores in results['summary']['celltype_scores'].items():
            print(f"  • {ct}: {scores['mean']:.1f} (range: {scores['min']}-{scores['max']})")
            
    except Exception as e:
        print(f"❌ Test Case 2 failed: {e}")
        import traceback
        traceback.print_exc()
    
    # Test Case 3: Custom models
    print("\n" + "="*60)
    print("📋 Test Case 3: Custom Model Configuration")
    print("="*60)
    
    try:
        # Use specific models
        custom_models = [
            "google/gemini-2.5-flash",
            "anthropic/claude-3.5-haiku"
        ]
        
        results = symphonyCompare(
            tissue="liver",
            celltypes=["Hepatocyte", "Kupffer cell"],
            marker_set="ALB, CYP3A4, CD163, MARCO",
            species="human",
            model_preset="custom",
            custom_models=custom_models,
            enable_discussion=False,  # No discussion for 2 models
            output_dir="./symphony_test_results",
            output_basename="liver_cells_custom",
            verbose=True
        )
        
        print(f"\n📊 Results Summary:")
        print(f"  • Models used: {len(custom_models)}")
        print(f"  • Consensus: {results['consensus'] or 'No consensus'}")
        
    except Exception as e:
        print(f"❌ Test Case 3 failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*70)
    print("✨ Symphony Compare demonstration completed!")
    print("="*70)


def show_feature_comparison():
    """Show comparison between old compareCelltypes and new symphonyCompare."""
    
    print("\n📊 Feature Comparison: compareCelltypes vs symphonyCompare")
    print("="*60)
    
    comparison = """
┌─────────────────────────┬──────────────────────┬──────────────────────┐
│ Feature                 │ compareCelltypes     │ symphonyCompare      │
├─────────────────────────┼──────────────────────┼──────────────────────┤
│ Function Name           │ compareCelltypes     │ symphonyCompare ✨   │
│ Multi-model Analysis    │ ✅ Yes               │ ✅ Yes               │
│ Parallel Processing     │ ✅ Yes               │ ✅ Yes               │
│ Auto Discussion Rounds  │ ✅ Yes (optional)    │ ✅ Yes (default on)  │
│ Consensus Detection     │ ✅ Basic             │ ✅ Advanced          │
│ Confidence Scores       │ ❌ No                │ ✅ Yes               │
│ Model Presets           │ ✅ 2 presets         │ ✅ 3+ presets        │
│ Custom Models           │ ✅ Yes               │ ✅ Yes               │
│ HTML Reports            │ ✅ Yes               │ ✅ Enhanced          │
│ CSV Output              │ ✅ Yes               │ ✅ Yes               │
│ Progress Tracking       │ ❌ Basic             │ ✅ Detailed          │
│ Summary Statistics      │ ❌ No                │ ✅ Yes               │
│ Output Organization     │ ❌ Basic             │ ✅ Advanced          │
│ API                     │ Multiple params      │ Unified interface    │
└─────────────────────────┴──────────────────────┴──────────────────────┘

Key Improvements in symphonyCompare:
• 🎼 Musical naming theme (Symphony, Movements, etc.)
• 📊 Returns confidence scores and summary statistics
• 🗂️ Better output file organization with timestamps
• 🎯 Configurable consensus threshold
• 📈 Detailed cell type score statistics (mean, std, range)
• 🎨 Enhanced progress messages and formatting
• 🔧 More intuitive parameter names
"""
    
    print(comparison)


if __name__ == "__main__":
    print("🎵 CASSIA Symphony Compare - Test Suite")
    print("="*50)
    
    # Show feature comparison
    show_feature_comparison()
    
    # Run the demonstration
    demo_symphony_compare()
    
    print("\n💡 Note: symphonyCompare is the enhanced version of compareCelltypes")
    print("   It provides a unified interface with automatic consensus building!")