COMPILER=gcc

LINK=../CLAPACK-3.2.1/lapack_MAC.a ../CLAPACK-3.2.1/blas_MAC.a ../CLAPACK-3.2.1/F2CLIBS/libf2c.a

LIB=-lc -lm


OPTIONS=-pg

pro_4comp: pro_4comp.c	
	$(COMPILER) $(LINK) $(LIB) pro_4comp.c -o pro_4comp.o 



pro_4comp_interaction: pro_4comp_interaction.c 	
	$(COMPILER) $(LINK) $(LIB) pro_4comp_interaction.c -o pro_4comp_interaction.o

pro_4comp_interaction_para: pro_4comp_interaction_para.c 	
	$(COMPILER) $(LINK) $(LIB) pro_4comp_interaction_para.c -o pro_4comp_interaction_para.o



pro_4comp_joint: pro_4comp_joint.c 	
	$(COMPILER) $(LINK) $(LIB) pro_4comp_joint.c -o pro_4comp_joint.o


pro_4comp_joint_para: pro_4comp_joint_para.c 	
	$(COMPILER) $(LINK) $(LIB) pro_4comp_joint_para.c -o pro_4comp_joint_para.o



pro_4comp_joint_para1: pro_4comp_joint_para1.c 	
	$(COMPILER) $(LINK) $(LIB) pro_4comp_joint_para1.c -o pro_4comp_joint_para1.o

