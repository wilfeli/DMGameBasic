DMGameBasic
===========

This code implements the basic Dynamic Macroeconomic (DM) Game model developed in the following paper:

Ekaterina Sinitskaya and Leigh Tesfatsion, "Macroeconomies as Constructively Rational Games," Journal of Economic Dynamics and Control, Vol. 61, 2015, 152-182.


Contacts:
kate.sinitskaya@gmail.com, tesfatsi@iastate.edu

Abstract:
Real-world decision-makers are forced to be locally constructive; that is, their decisions are necessarily constrained by their interaction networks, information, beliefs, and physical states. This study transforms an otherwise standard dynamic macroeconomic model into an open-ended dynamic game by requiring consumers and firms with intertemporal utility and profit objectives to be locally constructive. Tested locally-constructive decision processes for the consumers and firms range from simple reactive reinforcement learning to adaptive dynamic programming (ADP). Computational experiments are used to explore macroeconomic performance under alternative decision-process combinations relative to a social planner benchmark solution. A key finding is that simpler decision processes can outperform more sophisticated decision processes such as ADP. However, memory length permitting some degree of adaptive foresight is critical for good performance.

JEL Codes: B4, C6, C7, E03, E2

Keywords:  Macroeconomics; agent-based modeling; game theory; intertemporal optimization; learning; constructive rationality


Working paper preprint available at:

http://www2.econ.iastate.edu/tesfatsi/MacroConstructiveRationalityWP.SinitskayaTesfatsion.pdf


==============

To use:
---------

### Preferred way 

You would need Visual Studio 2015 for Windows 10
1. download solution file (and all other source and configuration code, cloning repository is enough for that)
2. compile to your preferred solution configuration (debug for ability to debug extensively, release if you want speed)
3. run


### Configuration
1. minimacro.exe requires two parameters: path to the configuration file, path to the save directory. Visual Studio solution/project file uses minimacro.ini file in the conf directory and saves results into Saves directory.
2. There are sample bash scripts in the conf\BashScripts directory that were used in the initial runs. They were utilized to generate .ini files for each configuration and run application file. They were tested on Mac (in 2014) and are included as a reference only. Do not expect them to work on Windows or Linux out of the box. 



### Manual compilation
1. download source code and compile using -std=c++11 flag (or similar)
2. You will need boost and Eigen headers (no install required) for compilation. Compatible headers are included with the source code, you might use them if you wish. 
3. minimacro.ini is a sample .ini file
