WannierBSE

🌐 Product Overview
WannierBSE is a high-performance computational suite engineered for solving the Bethe-Salpeter 
Equation (BSE) within the Wannier Tight-Binding (WTB) framework. Specifically optimized for 2D
materials and heterostructures, WannierBSE bridges the gap between first-principles accuracy 
and model-scale efficiency.

By leveraging MATLAB’s robust numerical engine, the package enables researchers to simulate 
excitonic spectra, dipoles, and light-matter interactions under complex dielectric environments.

🔗 Resources
Official Website: https://quantum.web.nycu.edu.tw/wannierbse

Source Code: https://github.com/WannierBSE/WannierBSE-MATLAB-

🛠 Technical Workflow & Integration
WannierBSE features a modular data pipeline designed to handle both raw generation and pre-computed
data ingestion.

1 Data Preparation:

Users can initiate calculations through two primary pathways:

Generation Route: Place wannier90_hr.dat in the /User_input/ folder. The solver will automatically
generate tight-binding bands, k-meshes, and dielectric functions, saving them in the /Precomputed_data/
folder for future use (e.g., v*_TB_WBSE.txt, c*_TB_WBSE.txt, kmesh_WBSE.txt, or epsilon_WBSE.txt).

Custom/Direct Route: Import user-defined files into /User_input/. The solver prioritizes files 
found here (e.g., v*_TB.txt, c*_TB.txt, kmesh.txt, or epsilon.txt) to bypass internal generation.

2 Configuration & Control:

The simulation environment is governed by dedicated control files located in the /Parameters/ directory:

   control.txt & structure.txt

   WF_centers.txt & kmesh_control.txt

   WTB_control.txt & dielectric_control.txt

3 Execution & Visualization:

Core Solver: Execute WBSE.m in MATLAB. The script utilizes an optimized, symmetry-aware Hamiltonian 
constructor and a high-performance parallelized solver.

Post-Processing: Use the auxiliary script EX_plot.m to visualize the exciton energy spectrum.

Data Management: All results—including matrices, energy spectra (Ex.mat), and wavefunctions 
(A.mat)—are automatically saved in the /Exciton_data/ folder.

📚 Citation

To support the continued development of WannierBSE, please cite the following references:

Primary Reference: Peng, G.-H. et al. "Distinctive Signatures of the Spin- and Momentum-Forbidden
Dark Exciton States in the Photoluminescence of Strained WSe2 Monolayers under Thermalization", 
Nano Lett. 19, 4, 2299–2312 (2019). https://doi.org/10.1021/acs.nanolett.8b04786

Non-local Dielectric Function: Li, W.-H. et al. "The Key Role of Non-Local Screening in the 
Environment-Insensitive Exciton Fine Structures of Transition-Metal Dichalcogenide Monolayers", 
Nanomaterials 13, 11, 1739 (2023). https://doi.org/10.3390/nano13111739

👥 Credits

Developed by the research group of Prof. Shun-Jen Cheng at National Yang Ming Chiao Tung 
University (NYCU).

Core Developers: Dr. Ping-Yuan Lo, Dr. Wei-Hua Li, Dr. Guan-Hao Peng, Dr. Jhen-Dong Lin, 
Dr. Vo Khuong Dien, Dr. Oscar Javier Gomez Sanchez, Mr. Ching-Hung Shih, Mr. Kun-Yi Lin, 
and Prof. Shun-Jen Cheng.

Copyright © 2026 The WannierBSE Development Team. All rights reserved.