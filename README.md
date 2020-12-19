# SpatialAllocationCOVID19


\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage[english]{babel}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{cite}
\usepackage{footnote}
\makesavenoteenv{tabular}
\usepackage{fancyref}
\usepackage{parskip}
\usepackage{url}
\usepackage{hyperref}

\usepackage[left=2cm,right=2cm,top=2cm,bottom=2cm]{geometry}

\title{Spatial Allocation of Scarce Vaccine and Antivirals for COVID-19 \\ \medskip \Large{ README File for Matlab/TOMLAB Code}}

\author{François M. Castonguay,\footnote{Department of Agricultural and Resource Economics, University of California, Davis, Davis, CA 95616, USA. Corresponding author. e-mail: \href{mailto: fcastonguay@ucdavis.edu}{fcastonguay@ucdavis.edu}}~ Julie C. Blackwood,\footnote{Department of Mathematics and Statistics, Williams College, Williamstown, MA 01267, USA.}~ Emily Howerton,\footnote{Department of Biology, Pennsylvania State University, University Park, PA 16802, USA.} \\ Katriona Shea,\footnotemark[3]~ Charles Sims,\footnote{Department of Economics, University of Tennessee, Knoxville, Knoxville, TN 37996, USA.}~ and James N. Sanchirico\footnote{Department of Environmental Science and Policy, University of California, Davis, Davis, CA 95616, USA and Resources for the Future, Washington DC 20036, USA.}}

\date{\nodate }

\begin{document}

\maketitle


\section{Problem Set-up}

\begin{itemize}
    \item``\texttt{covid\_Parameters\_months.m}"
    
    \begin{itemize}
        \item  This file contains all parameter values. It is called upon in ``Cost Effectiveness Main.m"; note that it is also copied in ``2. Problem Solving" so that the user can simply download this whole file to re-run the results.
    \end{itemize}
    
    \item ``\texttt{ConstraintOnDrug.m}"
    
    \begin{itemize}
        \item This is file is used to find the maximum level of infection of an uncontrolled outbreak that is used to determine the upper-bound constraint (i.e. 100\%) of the drug control. 
    \end{itemize}
    
    \item ``\texttt{InitialConditions.m}"
    
        \begin{itemize}
        \item This is the file used to find simulate out the outbreak to yield the initial conditions of our problem.
        \end{itemize}
    
    
\end{itemize}

\section{Problem Solving}

\begin{itemize}

    \item``\texttt{covid\_Parameters\_months.m}"
    
    \begin{itemize}
        \item  Same as above.
    \end{itemize}
    
    
    \item``\texttt{covid19\_Main.m}"
        \begin{itemize}
            \item This is the main file. The reader only needs to open this file to rerun our analysis. However, please make sure that ``\texttt{covid\_Parameters\_months.m}", ``\texttt{covid19\_Drugs.m}", and ``\texttt{covid19\_Vaccines.m}" are copied in the same folder. 
                \begin{itemize}
                    \item By choosing ``\texttt{omega=0}" and ``\texttt{ODE=1}", one will get the results from ``\texttt{workspacePN\_Gauss.mat}". Once yield, the workspace needs to be saved with the corresponding name.
                    \item By choosing ``\texttt{omega=0}" and ``\texttt{ODE=2}", one will get the results from ``\texttt{workspacePT\_Gauss.mat}". Once yield, the workspace needs to be saved with the corresponding name.
                    \item By choosing ``\texttt{omega=1/6}" and ``\texttt{ODE=1}", one will get the results from ``\texttt{workspaceST\_Gauss.mat}". Once yield, the workspace needs to be saved with the corresponding name.
                    \item By choosing  ``\texttt{omega=1/6}" and ``\texttt{ODE=2}", one will get the results from ``\texttt{workspaceSN\_Gauss.mat}". Once yield, the workspace needs to be saved with the corresponding name.
                \end{itemize}
            
        \end{itemize}
        

\pagebreak
        
    
    \item ``\texttt{covid19\_Drugs.m}"
            \begin{itemize}
            \item This file relates to the drug control of the problem. It details the state variables, the control variables, the ordinary differential equations describing the dynamics of the state variables, the constraints on state and controls variables across the different treatment scenarios (optimal and ad hoc), the objective of the central planning agency, how we generate our initial guess of the solution, and how we switch solvers if the chosen solver cannot find a solution.
            \end{itemize}
    
    
    \item ``\texttt{covid19\_Vaccines.m}"
            \begin{itemize}
            \item This file relates to the vaccine control of the problem. It details the state variables, the control variables, the ordinary differential equations describing the dynamics of the state variables, the constraints on state and controls variables across the different treatment scenarios (optimal and ad hoc), the objective of the central planning agency, how we generate our initial guess of the solution, and how we switch solvers if the chosen solver cannot find a solution.
            \end{itemize}
            
\end{itemize}

\section{Problem Analysis}

\begin{itemize}

    \item ``\texttt{ManuscriptFigures.m}"
    
        \begin{itemize}
        \item This file uses ``\texttt{workspacePN\_Gauss.mat}", ``\texttt{workspacePS\_Gauss.mat}", ``\texttt{workspaceSN\_Gauss.mat}", and ``\texttt{workspaceSS\_Gauss.mat}" to produce the figures in the manuscript.
        \end{itemize}
    
    \item ``\texttt{COVID\_Robustness\_Drugs.m}"
    
        \begin{itemize}
        \item This file uses ``\texttt{workspacePN\_Gauss.mat}", ``\texttt{workspacePS\_Gauss.mat}", ``\texttt{workspaceSN\_Gauss.mat}", and ``\texttt{workspaceSS\_Gauss.mat}" to produce the robustness checks with drugs.
        \end{itemize}
        
    \item ``\texttt{COVID\_Robustness\_Vaccines.m}"
    
        \begin{itemize}
        \item This file uses ``\texttt{workspacePN\_Gauss.mat}", ``\texttt{workspacePS\_Gauss.mat}", ``\texttt{workspaceSN\_Gauss.mat}", and ``\texttt{workspaceSS\_Gauss.mat}" to produce the robustness checks with vaccines.
        \end{itemize}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        \subsection{Workspace: No Pharmaceutical Intervention }


    \item ``\texttt{workspaceNoControls.mat}"
    
        \begin{itemize}
        \item This file contains the results for all cases where there are no pharmaceutical intervention undertaken.
        \end{itemize}


        
        \subsection{Workspaces: Main Results }


    \item ``\texttt{workspacePN\_Gauss.mat}"
    
        \begin{itemize}
        \item This file contains the main results when there is permanent immunity (P) and there is noncompliance to travel restrictions (N).
        \end{itemize}
    
    
    
    \item ``\texttt{workspacePT\_Gauss.mat}"
    
        \begin{itemize}
        \item This file contains the main results when there is permanent immunity (P) and there is compliance to travel restrictions (T).
        \end{itemize}
    
    
    \item ``\texttt{workspaceSN\_Gauss.mat}"
    
        \begin{itemize}
        \item This file contains the main results when there is six-month immunity (S) and there is noncompliance to travel restrictions (N).
        \end{itemize}
        

    \item ``\texttt{workspaceST\_Gauss.mat}"
    
        \begin{itemize}
        \item This file contains the main results when there is six-month immunity (S) and there is compliance to travel restrictions (T).
        \end{itemize}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                \subsection{Workspaces: Sensitivity Analysis of Vaccine Effectiveness }

            \item ``\texttt{workspace\_qV\_PN.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the vaccine effectiveness when there is permanent immunity (P) and there is noncompliance to travel restrictions (N).
        \end{itemize}
    
    
    
    \item ``\texttt{workspace\_qV\_PT.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the vaccine effectiveness when there is permanent immunity (P) and there is compliance to travel restrictions (T).
        \end{itemize}
    
    
    \item ``\texttt{workspace\_qV\_SN.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the vaccine effectiveness when there is six-month immunity (S) and there is noncompliance to travel restrictions (N).
        \end{itemize}
        

    \item ``\texttt{workspace\_qV\_ST.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the vaccine effectiveness when there is six-month immunity (S) and there is compliance to travel restrictions (T).
        \end{itemize}
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                            \subsection{Workspaces: Sensitivity Analysis of Drug Effectiveness }

            \item ``\texttt{workspace\_qD\_PN.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the drug effectiveness when there is permanent immunity (P) and there is noncompliance to travel restrictions (N).
        \end{itemize}
    
    
    
    \item ``\texttt{workspace\_qD\_PT.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the drug effectiveness when there is permanent immunity (P) and there is compliance to travel restrictions (T).
        \end{itemize}
    
    
    \item ``\texttt{workspace\_qD\_SN.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the drug effectiveness when there is six-month immunity (S) and there is noncompliance to travel restrictions (N).
        \end{itemize}
        

    \item ``\texttt{workspace\_qD\_ST.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the drug effectiveness when there is six-month immunity (S) and there is compliance to travel restrictions (T).
        \end{itemize}
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
 \subsection{Workspaces: Sensitivity Analysis of Workability Cost Parameter on Vaccines }

            
                        \item ``\texttt{workspace\_Cadj\_PN.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the workability cost parameter on vaccines when there is permanent immunity (P) and there is noncompliance to travel restrictions (N).
        \end{itemize}
    
    
    
    \item ``\texttt{workspace\_Cadj\_PT.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the workability cost parameter on vaccines when there is permanent immunity (P) and there is compliance to travel restrictions (T).
        \end{itemize}
    
    
    \item ``\texttt{workspace\_Cadj\_SN.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the workability cost parameter on vaccines when there is six-month immunity (S) and there is noncompliance to travel restrictions (N).
        \end{itemize}
        

    \item ``\texttt{workspace\_Cadj\_ST.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the workability cost parameter on vaccines when there is six-month immunity (S) and there is compliance to travel restrictions (T).
        \end{itemize}
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
     \subsection{Workspaces: Sensitivity Analysis of Workability Cost Parameter on Drugs }
        
            
                        \item ``\texttt{workspace\_Cadj\_PN\_Drug.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the workability cost parameter on drugs when there is permanent immunity (P) and there is noncompliance to travel restrictions (N).
        \end{itemize}
    
    
    
    \item ``\texttt{workspace\_Cadj\_PT\_Drug.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the workability cost parameter on drugs when there is permanent immunity (P) and there is compliance to travel restrictions (T).
        \end{itemize}
    
    
    \item ``\texttt{workspace\_Cadj\_SN\_Drug.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the workability cost parameter on drugs when there is six-month immunity (S) and there is noncompliance to travel restrictions (N).
        \end{itemize}
        

    \item ``\texttt{workspace\_Cadj\_ST\_Drug.mat}"
    
        \begin{itemize}
        \item This file contains the results of the sensitivity analysis of the workability cost parameter on drugs when there is six-month immunity (S) and there is compliance to travel restrictions (T).
        \end{itemize}
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




\end{itemize}




\end{document}


 
