DMGameBasic
===========

This code implements the basic Dynamic Macroeconomic (DM) Game model developed in the following working paper:
http://www2.econ.iastate.edu/tesfatsi/MacroConstructiveRationalityWP.SinitskayaTesfatsion.pdf


Ekaterina Sinitskaya and Leigh Tesfatsion, "Macroeconomies as Constructively Rational Games," Journal of Economic Dynamics and Control, Vol. 61, 2015, 152-182.


Contacts:
kate.sinitskaya@gmail.com, tesfatsi@iastate.edu

Abstract:
Real-world decision-makers are forced to be locally constructive; that is, their decisions are necessarily constrained by their interaction networks, information, beliefs, and physical states. This study transforms an otherwise standard dynamic macroeconomic model into an open-ended dynamic game by requiring consumers and firms with intertemporal utility and profit objectives to be locally constructive. Tested locally-constructive decision processes for the consumers and firms range from simple reactive reinforcement learning to adaptive dynamic programming (ADP). Computational experiments are used to explore macroeconomic performance under alternative decision-process combinations relative to a social planner benchmark solution. A key finding is that simpler decision processes can outperform more sophisticated decision processes such as ADP. However, memory length permitting some degree of adaptive foresight is critical for good performance.

JEL Codes: B4, C6, C7, E03, E2

Keywords:  Macroeconomics; agent-based modeling; game theory; intertemporal optimization; learning; constructive rationality

==============

To use:

1. download source code and compile using -std=c++11 flag (or similar)

2. You will need boost and Eigen headers (no install required) for compilation

3. minimacro.ini is a sample .ini file
