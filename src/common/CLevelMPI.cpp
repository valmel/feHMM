void CLevel::AssembleMaster() {
	int NumOfLeafEl = GetNumOfLeafEl();
	int par = 0;
	int tag = EMF_OK;
	double Center[DIM_OF_WORLD];
	double NextPrint;

	TRAVERSE_STACK *stack;
	static bool RestrictedElement;

	EL ** el_array= NULL;
	REAL * VolArray= NULL;
	int Proc;

	for (int i = 0; i < pAlberta->DimOfProb; i++)
		for (int j = 0; j < pAlberta->DimOfProb; j++) {
			clear_dof_matrix(pAlberta->Matrix[i][j]); //clear dof matrix in the beggining
			clear_dof_matrix(pAlberta->RestrictedMatrix[i][j]); //clear dof matrix in the beggining
		}

	el_array = (EL**)malloc((MPIInfo.Size+1)*sizeof(EL*));
	VolArray = (REAL*)malloc((MPIInfo.Size+1)*sizeof(REAL));

	NumOfLeafEl = GetNumOfLeafEl();
	double OneTenth = NumOfLeafEl/10.0;
	NextPrint = OneTenth; 

	stack = get_traverse_stack();

	pElInfo = traverse_first(stack, pAlberta->pMesh, -1, CALL_LEAF_EL
			| FILL_COORDS);
	for (Proc = 1; Proc < MPIInfo.Size && par < NumOfLeafEl; Proc++) { // start as many processes as possible
		// calculate center of element
		memset(Center, 0, pAlberta->DimOfWorld*sizeof(double)); // clear Center
		for (int i = 0; i < N_VERTICES(pAlberta->DimOfWorld); i++)
			for (int j = 0; j < pAlberta->DimOfWorld; j++)
				Center[j]+=pElInfo->coord[i][j];
		for (int i = 0; i < pAlberta->DimOfWorld; i++)
			Center[i]/=N_VERTICES(pAlberta->DimOfWorld);
		// save Volume and El
		VolArray[Proc] = el_det(pElInfo);
		el_array[Proc] = pElInfo->el;
		// move to next element
		pElInfo = traverse_next(stack, pElInfo);
		// send command to the slave
		tag = EMF_CALCULATE;
		MPI_Send(Center, pAlberta->DimOfWorld, MPI_DOUBLE, Proc, tag,
				MPI_COMM_WORLD);
		par = par + 1;
	}

	for (Proc = NumOfLeafEl + 1; Proc < MPIInfo.Size; Proc++) { // stop redundand processes (not probable)
		memset(Center, 0, pAlberta->DimOfWorld*sizeof(double));
		tag = EMF_STOP;
		MPI_Send(Center, pAlberta->DimOfWorld, MPI_DOUBLE, Proc, tag,
				MPI_COMM_WORLD);
	}

	int num_workers = (MPIInfo.Size-1) > NumOfLeafEl ? NumOfLeafEl
			: MPIInfo.Size-1;

	int pos;
	while (1) {
		MPI_Recv(pCommunicationMatrix, MatrixSize*MatrixSize, MPI_DOUBLE,
				MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &MPIInfo.Status);

		if (MPIInfo.Status.MPI_TAG == EMF_OK) {
			pos = 0;
			for (int i = 0; i < pAlberta->DimOfProb; i++)
				for (int j = 0; j < pAlberta->DimOfProb; j++)
					for (int k = 0; k < pAlberta->NumOfBasFct[i]; k++) {
						memcpy(pAlberta->LocalMatrix[i][j][k],
								pCommunicationMatrix + pos,
								pAlberta->NumOfBasFct[j]*sizeof(REAL));
						pos += pAlberta->NumOfBasFct[j];
					}

			RestrictedElement = true;
			for (int i = 0; i < pAlberta->DimOfProb; i++) {
				pAlberta->GetGlobDof[i](el_array[MPIInfo.Status.MPI_SOURCE],
						pAlberta->FEspace[i]->admin, pAlberta->Dofs[i]);
				// check wether it is interior element
				if (RestrictedElement)
					for (int j = 0; j < pAlberta->NumOfBasFct[i]; j++) {
						pAlberta->RestrictedDofs[i][j]
								=(DOF)pAlberta->pExtendedVertexIndex[i]->vec[pAlberta->Dofs[i][j]];
						if (pAlberta->RestrictedDofs[i][j] == -1) {
							RestrictedElement = false;
							break;
						}
					}
			}

			for (int i = 0; i<pAlberta->DimOfProb; i++)
				for (int j = 0; j<pAlberta->DimOfProb; j++) {
					if (pMatrixMask[i * pAlberta->DimOfProb + j]) {
						add_element_matrix(pAlberta->Matrix[i][j],
								VolArray[MPIInfo.Status.MPI_SOURCE],
								pAlberta->NumOfBasFct[i],
								pAlberta->NumOfBasFct[j], pAlberta->Dofs[i],
								pAlberta->Dofs[j],
								(const REAL**)pAlberta->LocalMatrix[i][j], nil);
						if (RestrictedElement)
							add_element_matrix(pAlberta->RestrictedMatrix[i][j],
									VolArray[MPIInfo.Status.MPI_SOURCE],
									pAlberta->NumOfBasFct[i],
									pAlberta->NumOfBasFct[j],
									pAlberta->RestrictedDofs[i],
									pAlberta->RestrictedDofs[j],
									(const REAL**)pAlberta->LocalMatrix[i][j],
									nil);
					}
				}
			//
			if (par < NumOfLeafEl) {
				// calculate center of element
				memset(Center, 0, pAlberta->DimOfWorld*sizeof(double)); // clear Center
				for (int i = 0; i < N_VERTICES(pAlberta->DimOfWorld); i++)
					for (int j = 0; j < pAlberta->DimOfWorld; j++)
						Center[j]+=pElInfo->coord[i][j];
				for (int i = 0; i < pAlberta->DimOfWorld; i++)
					Center[i]/=N_VERTICES(pAlberta->DimOfWorld);
				// save Volume and El
				VolArray[MPIInfo.Status.MPI_SOURCE] = el_det(pElInfo);
				el_array[MPIInfo.Status.MPI_SOURCE] = pElInfo->el;
				// move to next element
				pElInfo = traverse_next(stack, pElInfo);

				Proc = MPIInfo.Status.MPI_SOURCE;
				tag = EMF_CALCULATE;
				MPI_Send(Center, pAlberta->DimOfWorld, MPI_DOUBLE, Proc, tag,
						MPI_COMM_WORLD);
				par = par + 1;

				if (NextPrint <= par) {
					cout<<(int)10*NextPrint/OneTenth<<"%|";
					NextPrint += OneTenth;
					flush(cout);
				}

			} else {
				num_workers = num_workers - 1;
				Proc = MPIInfo.Status.MPI_SOURCE;
				memset(Center, 0, pAlberta->DimOfWorld*sizeof(double));
				tag = EMF_STOP;
				MPI_Send(Center, pAlberta->DimOfWorld, MPI_DOUBLE, Proc, tag,
						MPI_COMM_WORLD);

				if (num_workers == 0) {
					break;
				}
			}
		}
	}
	if (NextPrint <= par)
		cout<<"100%|";
	cout<<endl;
	pAlberta->AssembleProblemMatrix(pMatrixMask,
			(const DOF_MATRIX ***)pAlberta->Matrix);
	pAlberta->PostProcessMatrix(1.);
	if (el_array) {
		free(el_array);
		el_array = NULL;
	}
	if (VolArray) {
		free(VolArray);
		VolArray = NULL;
	}
	free_traverse_stack(stack);
}

void CLevel::AssembleSlave() {
	double Center[DIM_OF_WORLD];

	int tag, is, pos;
	TRAVERSE_STACK *stack;
	REAL bary[N_LAMBDA]= { 0 };

	stack = get_traverse_stack();

	while (1) {
		MPI_Recv(&Center, pAlberta->DimOfWorld, MPI_DOUBLE, 0, MPI_ANY_TAG,
				MPI_COMM_WORLD, &MPIInfo.Status);
		switch (MPIInfo.Status.MPI_TAG) {
		case EMF_OK:
			continue;
		case EMF_STOP:
			free_traverse_stack(stack);
			return;
		case EMF_CALCULATE:
			for (pElInfo = traverse_first(stack, pAlberta->pMesh, -1,
					CALL_LEAF_EL | FILL_COORDS); pElInfo; pElInfo
					= traverse_next(stack, pElInfo)) {
				is = world_to_coord(pElInfo, Center, bary);
				if (is == -1)
					break;
			}
			if (pElInfo == NULL) {
				printf("el_info not found\n");
				return;
			}
			for (int i = 0; i < pAlberta->DimOfProb; i++)
				pAlberta->GetGlobDof[i](pElInfo->el, pAlberta->FEspace[i]->admin,
						pAlberta->Dofs[i]);

			ComputeForElInfo(pElInfo);

			memset(pCommunicationMatrix, 0, MatrixSize*MatrixSize
					*sizeof(double));
			pos = 0;
			for (int i = 0; i < pAlberta->DimOfProb; i++)
				for (int j = 0; j < pAlberta->DimOfProb; j++)
					for (int k = 0; k < pAlberta->NumOfBasFct[i]; k++) {
						memcpy(pCommunicationMatrix + pos,
								pAlberta->LocalMatrix[i][j][k],
								pAlberta->NumOfBasFct[j]*sizeof(REAL));
						pos += pAlberta->NumOfBasFct[j];
					}
			tag = EMF_OK;
			MPI_Send(pCommunicationMatrix, MatrixSize*MatrixSize, MPI_DOUBLE,
					0, tag, MPI_COMM_WORLD);
			break;
		}
	}
}

void CLevel::SetUpRHSMaster() {
	int NumOfLeafEl = GetNumOfLeafEl();
	int par = 0;
	int tag = EMF_OK;
	double Center[DIM_OF_WORLD];
	double NextPrint;

	TRAVERSE_STACK *stack;

	EL **el_array = NULL;
	int Proc;

	el_array = (EL**)malloc((MPIInfo.Size+1)*sizeof(EL*));

	NumOfLeafEl = GetNumOfLeafEl();
	double OneTenth = NumOfLeafEl/10.0;
	NextPrint = OneTenth;

	for (int i = 0; i < pAlberta->DimOfProb; i++)
		dof_set(0.0, pAlberta->Fh[i]);

	stack = get_traverse_stack();
	pElInfo = traverse_first(stack, pAlberta->pMesh, -1, CALL_LEAF_EL
			| FILL_COORDS);
	for (Proc = 1; Proc < MPIInfo.Size && par < NumOfLeafEl; Proc++) { // start as many processes as possible
		// calculate center of element
		memset(Center, 0, pAlberta->DimOfWorld*sizeof(double)); // clear Center
		for (int i = 0; i < N_VERTICES(pAlberta->DimOfWorld); i++)
			for (int j = 0; j < pAlberta->DimOfWorld; j++)
				Center[j]+=pElInfo->coord[i][j];
		for (int i = 0; i < pAlberta->DimOfWorld; i++)
			Center[i]/=N_VERTICES(pAlberta->DimOfWorld);
		// save Volume and El
		el_array[Proc] = pElInfo->el;
		// move to next element
		pElInfo = traverse_next(stack, pElInfo);
		// send command to the slave
		tag = EMF_CALCULATE;
		MPI_Send(Center, pAlberta->DimOfWorld, MPI_DOUBLE, Proc, tag,
				MPI_COMM_WORLD);
		par = par + 1;
	}

	for (Proc = NumOfLeafEl + 1; Proc < MPIInfo.Size; Proc++) { // stop redundand processes (not probable)
		memset(Center, 0, pAlberta->DimOfWorld*sizeof(double));
		tag = EMF_STOP;
		MPI_Send(Center, pAlberta->DimOfWorld, MPI_DOUBLE, Proc, tag,
				MPI_COMM_WORLD);
	}

	int num_workers = (MPIInfo.Size-1) > NumOfLeafEl ? NumOfLeafEl
			: MPIInfo.Size-1;

	int pos;
	while (1) {
		MPI_Recv(pCommunicationVector, MatrixSize, MPI_DOUBLE, MPI_ANY_SOURCE,
				MPI_ANY_TAG, MPI_COMM_WORLD, &MPIInfo.Status);
		if (MPIInfo.Status.MPI_TAG == EMF_OK) {
			pos = 0;
			for (int i = 0; i<pAlberta->DimOfProb; i++) {
				memcpy(pAlberta->LocalFh[i], pCommunicationVector + pos,
						pAlberta->NumOfBasFct[i]*sizeof(REAL));
				pos += pAlberta->NumOfBasFct[i];
			}
			for (int i = 0; i < pAlberta->DimOfProb; i++)
				pAlberta->GetGlobDof[i](el_array[MPIInfo.Status.MPI_SOURCE],
						pAlberta->FEspace[i]->admin, pAlberta->Dofs[i]);
			for (int i = 0; i<pAlberta->DimOfProb; i++)
				for (int j = 0; j<pAlberta->NumOfBasFct[i]; j++) {
					pAlberta->Fh[i]->vec[pAlberta->Dofs[i][j]]
							+= pAlberta->LocalFh[i][j];
				}
			//
			if (par < NumOfLeafEl) {
				// calculate center of element
				memset(Center, 0, pAlberta->DimOfWorld*sizeof(double)); // clear Center
				for (int i = 0; i < N_VERTICES(pAlberta->DimOfWorld); i++)
					for (int j = 0; j < pAlberta->DimOfWorld; j++)
						Center[j]+=pElInfo->coord[i][j];
				for (int i = 0; i < pAlberta->DimOfWorld; i++)
					Center[i]/=N_VERTICES(pAlberta->DimOfWorld);
				// save Volume and El
				el_array[MPIInfo.Status.MPI_SOURCE] = pElInfo->el;
				// move to next element
				pElInfo = traverse_next(stack, pElInfo);

				Proc = MPIInfo.Status.MPI_SOURCE;
				tag = EMF_CALCULATE;
				MPI_Send(Center, pAlberta->DimOfWorld, MPI_DOUBLE, Proc, tag,
						MPI_COMM_WORLD);
				par = par + 1;

				if (NextPrint <= par) {
					cout<<(int)10*NextPrint/OneTenth<<"%|";
					NextPrint += OneTenth;
				}

			} else {
				num_workers = num_workers - 1;
				Proc = MPIInfo.Status.MPI_SOURCE;
				memset(Center, 0, pAlberta->DimOfWorld*sizeof(double));
				tag = EMF_STOP;
				MPI_Send(Center, pAlberta->DimOfWorld, MPI_DOUBLE, Proc, tag,
						MPI_COMM_WORLD);

				if (num_workers == 0) {
					break;
				}
			}
		}
	}
	if (NextPrint <= par)
		cout<<"100%|";
	printf("\n");
	if (el_array) {
		free(el_array);
		el_array = NULL;
	}
	int idx=0;
	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		for (int j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++) {
			pAlberta->ProblemRHS->vec[idx]=pAlberta->Fh[i]->vec[j];
			idx++;
		}
	}
	free_traverse_stack(stack);
}

void CLevel::SetUpRHSSlave() {
	double Center[DIM_OF_WORLD];

	int tag, is, pos;
	TRAVERSE_STACK *stack;
	REAL bary[N_LAMBDA]= { 0 };

	stack = get_traverse_stack();
	while (1) {
		MPI_Recv(&Center, pAlberta->DimOfWorld, MPI_DOUBLE, 0, MPI_ANY_TAG,
				MPI_COMM_WORLD, &MPIInfo.Status);
		switch (MPIInfo.Status.MPI_TAG) {
		case EMF_OK:
			continue;
		case EMF_STOP:
			free_traverse_stack(stack);
			return;
		case EMF_CALCULATE:
			for (pElInfo = traverse_first(stack, pAlberta->pMesh, -1,
					CALL_LEAF_EL | FILL_COORDS); pElInfo; pElInfo
					= traverse_next(stack, pElInfo)) {
				is = world_to_coord(pElInfo, Center, bary);
				if (is == -1)
					break;
			}
			if (pElInfo == NULL)
				exit(0);
			for (int i = 0; i < pAlberta->DimOfProb; i++)
				pAlberta->GetGlobDof[i](pElInfo->el, pAlberta->FEspace[i]->admin,
						pAlberta->Dofs[i]);
			SetUpRHSForOneElement(pElInfo);
			memset(pCommunicationVector, 0, MatrixSize*sizeof(double));
			pos = 0;
			for (int i = 0; i<pAlberta->DimOfProb; i++) {
				memcpy(pCommunicationVector + pos, pAlberta->LocalFh[i],
						pAlberta->NumOfBasFct[i]*sizeof(REAL));
				pos += pAlberta->NumOfBasFct[i];
			}
			tag = EMF_OK;
			MPI_Send(pCommunicationVector, MatrixSize, MPI_DOUBLE, 0, tag,
					MPI_COMM_WORLD);
			break;
		}
	}
}
