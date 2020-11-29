handle_read_draw_molecule_with_recentre ("A24t.pdb", 1)



delete_residue_range(0,'H',1,118)
delete_residue_range(0,'L',1,105)
renumber_residue_range(0,'L',106,300,2000)
renumber_residue_range(0,'H',119,306,3000)

change_chain_id(0,'H','L',1,3119,3300)
change_chain_id(0,'L','A',1,2106,3300)

newchain_info=merge_molecules([0,1], 0);

c=newchain_info[1][0];

change_chain_id(0,c,'A',1,335,526);

change_chain_id(0,c,'A',1,2001,2105);
change_chain_id(0,c,'A',1,3001,3118);

sort_chains(0);
sort_residues(0);
'''
