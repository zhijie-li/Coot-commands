set_nomenclature_errors_on_read ("autocorrect")
set_nomenclature_errors_on_read ("ignore")

set_refine_auto_range_step (2)
refine_residue_sphere_radius = 1.5
set_show_unit_cells_all (0)



def del_chains(imol,chains):
    if isinstance(chains,basestring):
        chains=list(chains)
    for c in chains:
        delete_chain(imol,c)
def sel(chain,res):
    return [[chain,i,''] for i in res]
    

def sel_range(chain,start,stop):
    return [[chain,i,''] for i in range(start,stop+1)]


handle_read_draw_molecule_with_recentre ("A.pdb", 1)
del_chains(0,'ABCDEF')
s1=sel_range('B',20,337)
s=[['B',5081,'']]+s1
refine_residues(0,s)

copy_residue_range(5,'L',12,'A',105,214)

delete_residue_range(0,'A',1,2000)

