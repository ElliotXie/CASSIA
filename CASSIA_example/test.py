#!/usr/bin/env python3
import sys
import os
import pandas as pd
import numpy as np
import time
import threading
import traceback
from pathlib import Path

openrouter_api_key = "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"
os.environ["OPENROUTER_API_KEY"] = openrouter_api_key


# Add the parent directory to sys.path to import CASSIA
sys.path.append(str(Path(__file__).parent.parent))

# Import CASSIA functions
from CASSIA_python.CASSIA.tools_function import runCASSIA_pipeline, runCASSIA_batch, set_api_key
from CASSIA_python.CASSIA.main_function_code import run_cell_type_analysis

# Model and provider configuration for all tests
MODEL = "deepseek/deepseek-chat-v3-0324"
PROVIDER = "openrouter"

# Set a timeout for operations (in seconds)
TIMEOUT = 1800  # 30 minutes (increased from 15 minutes)

def log_message(message, level="INFO"):
    """Log a message with timestamp and level"""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{level}] {message}")

class TimeoutError(Exception):
    """Exception raised when an operation times out"""
    pass

def create_sample_marker_data():
    """Create a sample marker dataset for testing with extremely long marker lists to trigger retries"""
    log_message("Creating sample marker data with extra long marker lists...")
    
    # Define some common cell types and their markers (with extremely extended marker lists)
    data = {
        'Cell Type': [
            'T Cell',
            'B Cell', 
            'Macrophage',
            'Neutrophil',
            'Endothelial Cell'
        ],
        'Markers': [
            'CD3D, CD3E, CD8A, CD4, TRAC, IL7R, CCR7, PTPRC, CD2, CD28, ICOS, CTLA4, PDCD1, CD27, LCK, ZAP70, LAT, CD247, IL2RA, CD69, SELL, S1PR1, CXCR3, CCR5, IFNG, IL2, TNF, GZMA, GZMB, PRF1, FOXP3, IL7, CD44, CD5, THEMIS, ITGAL, ITGB2, CD40LG, CD38, CD45RA, CD45RO, IL17A, IL17F, RORC, TBX21, GATA3, BCL6, STAT1, STAT3, STAT4, STAT5A, STAT5B, STAT6, IRF1, IRF4, ETS1, NFATC1, NFKB1, RELA, RUNX1, RUNX3, TOX, TOX2, LEF1, TCF7, ID2, ID3, MYB, EOMES, PRDM1, BCL2, BCL2L1, MCL1, CXCL13, CCL3, CCL4, CCL5, XCL1, XCL2, CXCL9, CXCL10, CXCL11, IFNA1, IFNA2, IFNB1, IL4, IL5, IL9, IL10, IL13, IL21, IL22, IL23A, IL12A, IL12B, IL15, IL18, IL1A, IL1B, IL6, CD8B, CD4, CD3G, CD3D, CD3E, TRDC, TRGC1, TRGC2, TRBC1, TRBC2, TCRA, TCRB, TRAV, TRBV, TRDV, TRGV, TRAJ, TRBJ, TRDJ, TRGJ, CSF2, TNFSF11, FASLG, CD70',
            'CD19, MS4A1, CD79A, CD79B, IGHM, IGHD, CD27, CD20, PAX5, EBF1, VPREB1, IGLL1, BLNK, BTK, CD22, CD24, FCER2, CR2, CD40, AICDA, RAG1, RAG2, IGKC, IGLC1, IGHG1, IGHG2, IGHG3, IGHG4, IGHA1, IGHA2, IGHE, MZB1, XBP1, PRDM1, IRF4, CXCR5, CCR6, BAFF, APRIL, TNFRSF13B, TNFRSF13C, TNFRSF17, SDC1, CD5, CD23, CD38, CD138, BLIMP1, BCL6, BACH2, MYC, IRF8, SYK, PTPRC, CD45RA, CD45RO, CIITA, HLA-DRA, HLA-DRB1, HLA-DQA1, HLA-DQB1, CXCR4, CXCL12, CCL19, CCL21, LYN, PIK3CD, PIK3R1, AKT1, AKT2, AKT3, FOXO1, FOXO3, FOXP1, BCL2, BCL2L1, MCL1, BCL2A1, BAX, BAK1, BIM, BAD, BID, FCGR2B, FCRL1, FCRL2, FCRL3, FCRL4, FCRL5, FCRLA, CD72, CD81, CD82, CD83, BANK1, BLK, GRB2, PLCG2, VAV1, CARD11, MALT1, BCL10, REL, RELA, NFKB1, NFKB2, TRAF1, TRAF2, TRAF3, MAP3K7, MAPK1, MAPK3, MAPK8, MAPK14, JAK1, JAK2, JAK3, TYK2, STAT1, STAT3, STAT5A, STAT5B, STAT6, IFI16, TLR7, TLR9, CD180, IL4R, IL13RA1, IL21R, CD19, CR2, CD27, CD40, CD79A, CD79B, IGHM, IGHD, IGHG, IGHA, IGHE, SHM, CSR, SLC7A7',
            'CD68, CD14, CSF1R, CD163, FCGR3A, MARCO, MSR1, ITGAM, ITGAX, ADGRE1, CD86, CD80, MRC1, MSR1, TREM2, CLEC7A, CD209, SIGLEC1, VSIG4, MERTK, STAB1, C1QA, C1QB, C1QC, TGFB1, PDGFB, VEGFA, IL10, IL1B, IL6, TNF, CCL2, CCL3, CCL4, CXCL8, CXCL9, CXCL10, HLA-DRA, HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DPA1, HLA-DPB1, ITGAL, ITGB2, FCGR1A, FCGR2A, FCGR2B, FCGR2C, CD55, CD59, CD44, CD74, CTSS, CTSL, CTSB, CTSD, CTSZ, LYZ, AIF1, APOE, TYROBP, FCER1G, NR1H3, PPARG, PPARA, PPARD, ABCA1, ABCG1, SCARB1, LPL, TGFBR1, TGFBR2, SMAD2, SMAD3, SMAD4, SMAD7, IFNGR1, IFNGR2, STAT1, IRF1, IRF3, IRF5, IRF7, IRF8, TLR1, TLR2, TLR3, TLR4, TLR5, TLR6, TLR7, TLR8, TLR9, MYD88, TIRAP, TRIF, TRAM, NOD1, NOD2, NLRP3, NLRC4, CASP1, IL1A, IL1B, IL18, IL33, S100A8, S100A9, S100A12, CD274, PDCD1LG2, TNFSF10, TNFSF9, TNFSF14, CIITA, SOCS1, SOCS3, CX3CR1, CCR2, CCR5, CXCR3, CXCR4, TGFBR, MAFB, ARG1, NOS2, CD38, CD36, CD9, CD81',
            'ELANE, MPO, FCGR3B, CSF3R, S100A8, S100A9, CXCR2, CXCR1, ITGAM, ITGB2, FPR1, FPR2, CEACAM8, CEACAM3, CD177, OLFM4, LTF, LCN2, BPI, DEFA1, DEFA3, DEFA4, CAMP, CRISP3, AZU1, PRTN3, CTSG, MMP8, MMP9, CXCL1, CXCL2, CXCL8, FCGR2A, FCGR2B, CD66B, CD11B, CD18, CLEC7A, TLR2, TLR4, NOD2, NLRP3, IL1B, TNF, NCF1, NCF2, NCF4, CYBB, CYBA, IL1R1, IL1R2, IL1RN, IL1RAP, IL1RL1, IL1RL2, IL18R1, IL18RAP, TNFRSF1A, TNFRSF1B, CD55, CD59, CR1, ITGAL, SLC11A1, FCGR1A, FCGR1B, FCGR3A, PGLYRP1, PGLYRP2, PGLYRP3, PGLYRP4, NLRC4, PYCARD, CASP1, CASP4, CASP5, IL1A, IL6, CCL3, CCL4, CCL2, CXCL1, CXCL2, CXCL3, CXCL5, CXCL6, CXCL7, CXCL9, CXCL10, CXCL11, CCR1, CCR2, CCR3, CXCR1, CXCR2, S100A12, ALOX5, ALOX5AP, LTA4H, PTGS1, PTGS2, ANXA1, ANXA3, CD14, CD15, IFNGR1, IFNGR2, CD64, CD16, CD32, CD35, CD88, CD11C, ICAM1, SELPLG, CD44, PECAM1, SELL, SIGLEC5, SIGLEC9, CD33, LAMP1, LAMP2, BCL2, BCL2L1, MCL1, BAX, CASP3, TRADD, TRAF2, RIPK1',
            'PECAM1, CDH5, VWF, CLDN5, KDR, TEK, ICAM1, VCAM1, SELE, SELP, NOS3, VEGFR2, TIE1, TIE2, ANGPT1, ANGPT2, VEGFA, VEGFB, VEGFC, VEGFD, FLT1, NRP1, NRP2, ROBO4, DLL4, JAG1, JAG2, NOTCH1, NOTCH4, HEY1, HEY2, EPHB4, EFNB2, GJA4, GJA5, CXCL12, CXCR4, EDN1, EDNRB, ACTA2, PDGFRB, CNN1, TAGLN, MYH11, ABCG2, LYVE1, PROX1, FLT4, PDPN, ITGA9, CCL21, ACKR3, CCL2, CCL5, CX3CL1, CD34, CD31, ICAM2, MCAM, CDH2, PODXL, THBD, THBS1, PROCR, F3, VWF, CD36, PLVAP, CAV1, AQP1, ADGRF5, ESAM, CD44, CD47, CD93, CD109, CD146, CD105, ENG, ACVRL1, SOX17, SOX18, ERG, FLI1, TAL1, GATA2, FOXP1, FOXF1, FOXC1, FOXC2, FOXO1, TWIST1, SNAI1, SNAI2, ZEB1, ZEB2, TGFB1, TGFB2, TGFB3, TGFBR1, TGFBR2, BMPR2, ACVR2A, SMAD1, SMAD2, SMAD3, SMAD4, SMAD5, BMP2, BMP4, BMP6, BMP9, ID1, ID2, ID3, HIF1A, EPAS1, CTGF, CCN1, CCN2, CCN3, MMP1, MMP2, MMP9, MMP14, TIMP1, TIMP2, SERPINE1, SERPINA1, SERPING1, ACE, ACE2, SOD1, CAT, GPX1, NOX1, NOX4, DUOX1, DUOX2, NQO1, NFE2L2, HMOX1'
        ]
    }
    
    # Convert to DataFrame
    df = pd.DataFrame(data)
    
    # Save to CSV
    csv_path = Path(__file__).parent / 'test_markers.csv'
    df.to_csv(csv_path, index=False)
    log_message(f"Sample marker data with extra long lists saved to {csv_path}")
    
    # Return the path as a string
    return str(csv_path.absolute())

class ProgressMonitor(threading.Thread):
    """Thread to monitor progress and log status periodically"""
    def __init__(self, interval=30):
        super().__init__()
        self.interval = interval
        self.daemon = True
        self.stopped = threading.Event()
        self.start_time = time.time()
        
    def run(self):
        while not self.stopped.wait(self.interval):
            elapsed = time.time() - self.start_time
            log_message(f"Still running... (elapsed time: {elapsed:.1f} seconds)", "DEBUG")
            
    def stop(self):
        self.stopped.set()

class TimeoutMonitor(threading.Thread):
    """Thread to monitor for timeout and raise an exception"""
    def __init__(self, timeout_seconds, target_thread):
        super().__init__()
        self.daemon = True
        self.timeout_seconds = timeout_seconds
        self.target_thread = target_thread
        self.stopped = threading.Event()

    def run(self):
        if not self.stopped.wait(self.timeout_seconds):
            if self.target_thread.is_alive():
                log_message(f"Operation timed out after {self.timeout_seconds} seconds", "ERROR")
                # In Python we can't force-kill a thread safely, so we'll
                # set a flag to check later
                self.target_thread._timed_out = True

    def stop(self):
        self.stopped.set()

def run_with_timeout(func, timeout_seconds, *args, **kwargs):
    """Run a function with a timeout"""
    result = [None]
    error = [None]
    finished = [False]
    
    def wrapped_func():
        try:
            result[0] = func(*args, **kwargs)
        except Exception as e:
            error[0] = e
        finally:
            finished[0] = True
    
    # Create and start the worker thread
    worker = threading.Thread(target=wrapped_func)
    worker._timed_out = False
    worker.daemon = True
    worker.start()
    
    # Create and start the timeout monitor
    timeout_monitor = TimeoutMonitor(timeout_seconds, worker)
    timeout_monitor.start()
    
    # Wait for the worker to finish
    worker.join(timeout_seconds + 1)  # Give a small buffer
    timeout_monitor.stop()
    
    # Check if we timed out
    if worker._timed_out or not finished[0]:
        raise TimeoutError(f"Operation timed out after {timeout_seconds} seconds")
    
    # Check if there was an error
    if error[0] is not None:
        raise error[0]
    
    return result[0]

def test_cassia_batch_with_retry():
    """Test the CASSIA batch functionality with retry mechanism"""
    
    # Create test data
    marker_path = create_sample_marker_data()
    
    log_message("\n=== Testing CASSIA batch with retry mechanism ===")
    log_message(f"Using model: {MODEL}")
    log_message(f"Using provider: {PROVIDER}")
    log_message(f"Using marker file: {marker_path}")
    
    # Verify the marker file exists
    if not os.path.exists(marker_path):
        log_message(f"ERROR: Marker file not found at {marker_path}", "ERROR")
        return
    
    # Make sure your API keys are set
    if not os.environ.get("OPENROUTER_API_KEY"):
        log_message("WARNING: OPENROUTER_API_KEY not found in environment.", "WARNING")
        log_message("Please set OPENROUTER_API_KEY or uncomment and use the set_api_key() function in the script.", "WARNING")
        return
    
    # Start progress monitor
    monitor = ProgressMonitor()
    monitor.start()
    
    start_time = time.time()
    
    # Run batch analysis with retry mechanism
    try:
        log_message("Starting CASSIA batch processing...")
        
        def run_cassia_batch():
            return runCASSIA_batch(
                marker=marker_path,
                output_name="test_results2",
                tissue="lung",
                species="human",
                additional_info="Sample test data with immune cells",
                max_workers=5,  # Use fewer workers for testing
                max_retries=1,  # Increased to 3 retries for better testing
                model=MODEL,
                provider=PROVIDER
            )
        
        # Run with timeout
        results = run_with_timeout(run_cassia_batch, TIMEOUT)
        
        elapsed = time.time() - start_time
        log_message(f"\n✅ CASSIA batch completed successfully in {elapsed:.1f} seconds!")
        log_message(f"Results saved to test_results_full.csv and test_results_summary.csv")
        
    except TimeoutError:
        log_message(f"\n⚠️ Operation timed out after {TIMEOUT} seconds", "ERROR")
    except Exception as e:
        elapsed = time.time() - start_time
        log_message(f"\n❌ Error running CASSIA batch after {elapsed:.1f} seconds: {str(e)}", "ERROR")
        log_message(f"Traceback: {traceback.format_exc()}", "ERROR")
    finally:
        # Stop monitor
        monitor.stop()

def test_cassia_pipeline():
    """Test the full CASSIA pipeline with retry functionality"""
    
    # Create test data
    marker_path = create_sample_marker_data()
    
    log_message("\n=== Testing full CASSIA pipeline with retry mechanism ===")
    log_message(f"Using model: {MODEL}")
    log_message(f"Using provider: {PROVIDER}")
    log_message(f"Using marker file: {marker_path}")
    
    # Verify the marker file exists
    if not os.path.exists(marker_path):
        log_message(f"ERROR: Marker file not found at {marker_path}", "ERROR")
        return
    
    # Ensure API keys are set in the environment
    if not os.environ.get("OPENROUTER_API_KEY"):
        log_message("WARNING: OPENROUTER_API_KEY not found in environment.", "WARNING")
        log_message("Please set OPENROUTER_API_KEY or uncomment and use the set_api_key() function in the script.", "WARNING")
        return
    
    # Start progress monitor
    monitor = ProgressMonitor()
    monitor.start()
    
    start_time = time.time()
    
    # Run the pipeline
    try:
        log_message("Starting CASSIA pipeline...")
        
        def run_cassia_pipeline():
            return runCASSIA_pipeline(
                output_file_name="pipeline_test2",
                tissue="lung",
                species="human",
                marker_path=marker_path,
                max_workers=5,  # Use fewer workers for testing
                max_retries=1,  # Set to 1 retry
                annotation_model=MODEL,
                annotation_provider=PROVIDER,
                score_model=MODEL,
                score_provider=PROVIDER,
                annotationboost_model=MODEL,
                annotationboost_provider=PROVIDER,
                additional_info="Sample test data with immune cells"
            )
        
        # Run with timeout
        run_with_timeout(run_cassia_pipeline, TIMEOUT)
        
        elapsed = time.time() - start_time
        log_message(f"\n✅ CASSIA pipeline completed successfully in {elapsed:.1f} seconds!")
        
    except TimeoutError:
        log_message(f"\n⚠️ Operation timed out after {TIMEOUT} seconds", "ERROR")
    except Exception as e:
        elapsed = time.time() - start_time
        log_message(f"\n❌ Error running CASSIA pipeline after {elapsed:.1f} seconds: {str(e)}", "ERROR")
        log_message(f"Traceback: {traceback.format_exc()}", "ERROR")
    finally:
        # Stop monitor
        monitor.stop()

if __name__ == "__main__":
    log_message("CASSIA Test Script")
    log_message("==================")
    log_message(f"Testing with model: {MODEL}")
    log_message(f"Testing with provider: {PROVIDER}")
    
    # Choose which test to run
    #test_cassia_batch_with_retry()
    
    # Run the full pipeline test
    test_cassia_pipeline()
    
    log_message("\nTest completed!")
