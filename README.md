DMGameBasic
===========

This code implements the basic Dynamic Macroeconomic (DM) Game model developed in the following working paper:
http://www.econ.iastate.edu/tesfatsi/MacroConstructiveRationalityWP.SinitskayaTesfatsion.pdf

Title: 

"Macroeconomies as Constructively Rational Games"

Authors: 

Ekaterina Sinitskaya and Leigh Tesfatsion 
Department of Economics, Iowa State University 
Ames, Iowa 50011-1070
kate.sinitskaya@gmail.com, tesfatsi@iastate.edu

Abstract:

Real-world decision-makers are forced to be locally constructive, in the sense that their actions are constrained by the interaction networks, limited information, and computational capabilities at their disposal.  This study poses the following question: Suppose utility-seeking consumers and profit-seeking firms in an otherwise standard dynamic macroeconomic model are required to be locally constructive decision-makers, unaided by the external imposition of global coordination conditions.  What combinations of locally constructive decision rules result in good macroeconomic performance relative to a social planner benchmark model, and what are the game-theoretic properties of these decision-rule combinations?  We begin our investigation of this question by specifying locally constructive decision rules for the consumers and firms that range from simple fixed behaviors to sophisticated approximate dynamic programming algorithms.  We then use computational experiments to explore macroeconomic performance under alternative decision-rule combinations.  A key finding is that simpler rules can outperform more sophisticated rules, but that forward-looking behavior coupled with a relatively long memory permitting past observations to inform current decision-making is critical for good performance.

JEL Codes: B4, C6, C7, E03, E2

Keywords:  Macroeconomics, agent-based, game, learning, stochastic optimization

==============
To use:

1. download source code and compile using -std=c++11 flag (or similar)

2. You will need boost and Eigen headers (no install required) for compilation

3. minimacro.ini is a sample .ini file
