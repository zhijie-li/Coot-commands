

set_nomenclature_errors_on_read ("autocorrect")
set_nomenclature_errors_on_read ("ignore")


set_refine_auto_range_step (2)
set_font_size (3)
#set_graphics_window_size (1920, 1000)
#set_graphics_window_position (0, 0)
refine_residue_sphere_radius = 1.5

#set_model_fit_refine_dialog_position(800, 100)
set_display_control_dialog_position (700, 300)
#set_go_to_atom_window_position (1, 78)


set_show_unit_cells_all (0)



#add "Glyco" menu
#add_module_carbohydrate_gui()

add_key_binding ("Zalman Stereo", "S", lambda:  zalman_stereo_mode())
add_key_binding ("Undo Zalman Stereo", "M", lambda: mono_mode())

add_key_binding("Refine Active Residue", "r", lambda: manual_refine_residues(0))
add_key_binding("Refine Active Residue AA", "x", lambda: refine_active_residue())

add_key_binding("Autofit Rotamer", "j", lambda: auto_fit_rotamer_active_residue())

add_key_binding("Pepflip", "q", lambda: pepflip_active_residue())
#add_key_binding("Eigen-flip Ligand", "e", lambda: flip_active_ligand())



add_key_binding ("Refine Auto Range", "A", lambda:  key_binding_func_Z1())
add_key_binding("Refine residue in a sphere", "R",lambda: sphere_refine(refine_residue_sphere_radius))
add_key_binding("Neighbours Refine", "H", lambda: key_binding_func_21())
add_key_binding("Jiggle Fit", "J", lambda: key_binding_func_7())
add_key_binding("Triple Refine", "T", lambda: manual_refine_residues(1))



add_key_binding("Go To Blob", "g", lambda: blob_under_pointer_to_screen_centre())

add_key_binding("Ca + Ligands", "B", lambda:   graphics_to_ca_plus_ligands_representation (0))
add_key_binding("Hydrogens off", "[", lambda: set_draw_hydrogens(0, 0))
add_key_binding("Hydrogens on", "]", lambda: set_draw_hydrogens(0, 1))

add_key_binding("Regularize Residues", "b", lambda: key_binding_func_6())
add_key_binding("Delete this water", "D", lambda: delete_atom(*active_residue()))
add_key_binding ("Add Water", "w", lambda: place_typed_atom_at_pointer ("Water"))
add_key_binding("Add terminal residue", "y", lambda: key_binding_func_1())
add_key_binding("Fill Partial", "k", lambda: key_binding_func_2())


def match_CA(s,e,ch,imolref,imolmove):
    simple_lsq_match(s,e,ch,imolref,s,e,ch,imolref, 'CA')
def match_main(s,e,ch,imolref,imolmove):
    simple_lsq_match(s,e,ch,imolref,s,e,ch,imolref, 'main')
def match_all(s,e,ch,imolref,imolmove):
    simple_lsq_match(s,e,ch,imolref,s,e,ch,imolref, 'all')

def hydrogen_switch(tf):
    for imol in range(20):        set_draw_hydrogens(imol,0)

def del_chains(imol,chains):
    if isinstance(chains,basestring):
        chains=list(chains)
    for c in chains:
        delete_chain(imol,c)

#def add_linked_residue(imol, chain_id, resno, ins_code, new_residue_comp_id, link_type, n_trials):
#add_nag(0,'A',17,5001)

Nglyc=[[0,"NAG", "NAG-ASN"],[1,"NAG", "BETA1-4"],[2,"BMA", "BETA1-4"],[3,"MAN", "ALPHA1-3"],[3,"MAN", "ALPHA1-6"]]

'''
def add_nag(aa_imol, aa_chain_id, aa_res_no,serial,aa_ins_code=''):
    delete_residue_range(aa_imol, aa_chain_id,serial,serial)
#    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
    newnag=add_linked_residue(aa_imol, aa_chain_id, aa_res_no, aa_ins_code, 'NAG', 'NAG-ASN',1)
    #[True, 'A', 5183, '']
    n=newnag[2]
    if newnag[0]:
        renumber_residue_range(aa_imol, aa_chain_id,n,n,serial-n)
        newnag[2]=serial
    return(newnag)
'''
def man3_here(res2):
    with UsingActiveAtom() as [i,ch,res1, aa_ins_code, aa_atom_name, aa_alt_conf]:
        add_Man3(i, ch, res1,res2)

def asnnag_here(res2):
    with UsingActiveAtom() as [i,ch,res1, aa_ins_code, aa_atom_name, aa_alt_conf]:
        add_Nglyc_monomer(i, ch, res1,res2,"NAG", "NAG-ASN")

def nag14_here():
    with UsingActiveAtom() as [i,ch,res1, aa_ins_code, aa_atom_name, aa_alt_conf]:
        add_Nglyc_monomer(i, ch, res1,res1+1,"NAG", "BETA1-4")

def bman_here():
    with UsingActiveAtom() as [i,ch,res1, aa_ins_code, aa_atom_name, aa_alt_conf]:
        add_Nglyc_monomer(i, ch, res1,res1+1,"BMA", "BETA1-4")

def add_Man3(i, ch, res1,res2,aa_ins_code=''):
    add_Nglyc_monomer(i, ch, res1,res2,"NAG", "NAG-ASN")
    add_Nglyc_monomer(i, ch, res2,res2+1,"NAG", "BETA1-4")
    add_Nglyc_monomer(i, ch, res2+1,res2+2,"BMA", "BETA1-4")


def add_2nag(i, aa_chain_id, aa_res_no,serial,aa_ins_code=''):
    add_Nglyc_monomer(i, ch, res1,res2,"NAG", "NAG-ASN")
    add_Nglyc_monomer(i, ch, res2,res2+1,"NAG", "BETA1-4")
#    delete_residue_range(aa_imol, aa_chain_id,serial,serial+1)
#    add_nag(aa_imol, aa_chain_id, aa_res_no,serial)
#    add_14nag(aa_imol, aa_chain_id, serial,serial+1)

def add_asnnag(i, ch,res1,res2):
    add_Nglyc_monomer(i, ch,res1,res2, "NAG", "NAG-ASN")
'''
def add_14nag(aa_imol, aa_chain_id, aa_res_no,serial):
    delete_residue_range(aa_imol, aa_chain_id,serial,serial)
    newnag=add_linked_residue(aa_imol, aa_chain_id, aa_res_no, '', 'NAG', "BETA1-4",1)
    n=newnag[2]
    if newnag[0]:
        renumber_residue_range(aa_imol, aa_chain_id,n,n,serial-n)
        newnag[2]=serial
    return newnag

'''
def add_Nglyc_monomer(aa_imol, aa_chain_id, res1,res2,sugar,link):
    delete_residue_range(aa_imol, aa_chain_id,res2,res2)
    newnag=add_linked_residue(aa_imol, aa_chain_id, res1, '', sugar,link,1)
    n=newnag[2]
    if newnag[0]:
        b=renumber_residue_range(aa_imol, aa_chain_id,n,n,res2-n)
        print(b)
        newnag[2]=res2
    return newnag

'''
def add_bman(aa_imol, aa_chain_id, aa_res_no,serial):
    delete_residue_range(aa_imol, aa_chain_id,serial,serial)
    newnag=add_linked_residue(aa_imol, aa_chain_id, aa_res_no, '', 'BMA', "BETA1-4",1)
    n=newnag[2]
    if newnag[0]:
        renumber_residue_range(aa_imol, aa_chain_id,n,n,serial-n)
        newnag[2]=serial
    return newnag
'''
def key_binding_func_Z1():
    active_atom = active_residue()
    if (not active_atom):
        add_status_bar_text("No active residue")
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        refine_auto_range(imol, chain_id,
												res_no,alt_conf)


def key_binding_func_7():
    active_atom = active_residue()
    if (not active_atom):
        add_status_bar_text("No active residue")
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        fit_to_map_by_random_jiggle(imol, chain_id,
												res_no,ins_code,100,1.0)

def key_binding_func_6():
    active_atom = active_residue()
    if (not active_atom):
        add_status_bar_text("No active residue")
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        regularize_zone(imol, chain_id,
                        res_no - 1, res_no + 1,
                        alt_conf)





def key_binding_func_1():
    active_atom = active_residue()
    if (not active_atom):
        print "No active atom"
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        add_terminal_residue(imol, chain_id, res_no, "auto", 1)


def key_binding_func_2():
    active_atom = active_residue()
    if (not active_atom):
        print "No active atom"
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        fill_partial_residue(imol, chain_id, res_no, ins_code)






def key_binding_func_21():
    if not valid_map_molecule_qm(imol_refinement_map()):
        info_dialog("Must set the refinement map")
    else:
        # not using active atom
        active_atom = active_residue()
        if (not active_atom):
            add_status_bar_text("No active residue")
        else:
            imol      = active_atom[0]
            chain_id  = active_atom[1]
            res_no    = active_atom[2]
            ins_code  = active_atom[3]
            atom_name = active_atom[4]
            alt_conf  = active_atom[5]

            rc_spec = [chain_id, res_no, ins_code]
            ls = residues_near_residue(imol, rc_spec, 4)
            with_auto_accept([refine_residues, imol, [rc_spec] + ls])


def stepped_refine_protein_sphere(imol, res_step = 2, rad =refine_residue_sphere_radius):

    import types
    from types import IntType

    imol_map = imol_refinement_map()
    if (not valid_map_molecule_qm(imol_map)):
        info_dialog("Oops, must set map to refine to")
    else:
        def refine_func(chain_id, res_no):
            #print "centering on ",chain_id,res_no," CA"
            set_go_to_atom_chain_residue_atom_name(chain_id, res_no, "CA")
            #rotate_y_scene(30, 0.3) # n-frames frame-interval(degrees)
            #refine_auto_range(imol, chain_id, res_no, "")
            sphere_refine(radius=rad, expand=False)
            accept_regularizement()
            #rotate_y_scene(30,0.3)
        stepped_refine_protein_with_refine_func_sphere(imol, refine_func, res_step)
def stepped_refine_protein_with_refine_func_sphere(imol, refine_func, res_step):

    set_go_to_atom_molecule(imol)
    make_backup(imol)
    backup_mode = backup_state(imol)
    alt_conf = ""
    imol_map = imol_refinement_map()
    replacement_state = refinement_immediate_replacement_state()

    if (imol_map == -1):
        # actually shouldnt happen as we set check the map earlier...
        add_status_bar_text("Oops.  Must set a map to fit")
    else:
        turn_off_backup(imol)
        set_refinement_immediate_replacement(1)
        res_step = int(res_step)
        if (res_step <= 1):
            set_refine_auto_range_step(1)
            res_step = 1
        else:
            set_refine_auto_range_step(int(res_step / 2))

        for chain_id in chain_ids(imol):
            n_residues = chain_n_residues(chain_id,imol)
            print "There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id}

            for serial_number in range(0, n_residues, res_step):
                res_name = resname_from_serial_number(imol, chain_id, serial_number)
                res_no = seqnum_from_serial_number(imol, chain_id, serial_number)
                ins_code = insertion_code_from_serial_number(imol, chain_id, serial_number)
                if ins_code is not None:
                    print "centering on ", chain_id, res_no, " CA"
                    refine_func(chain_id, res_no)

        if (replacement_state == 0):
            set_refinement_immediate_replacement(0)
        if (backup_mode == 1):
            turn_on_backup(imol)



def stepped_sphere_refine_function(res_spec, imol_map, use_rama = False):

    imol     = res_spec[0]
    chain_id = res_spec[1]
    res_no   = res_spec[2]
    ins_code = res_spec[3]
    current_steps_per_frame = dragged_refinement_steps_per_frame()
    current_rama_state = refine_ramachandran_angles_state()
    set_dragged_refinement_steps_per_frame(400)

    for alt_conf in residue_alt_confs(imol, chain_id, res_no, ins_code):
        if use_rama:
            set_refine_ramachandran_angles(1)
        res_name = residue_name(imol, chain_id, res_no, ins_code)
        if (not res_name == "HOH"):
#            print "centering on", chain_id, res_no, "CA"
            set_go_to_atom_chain_residue_atom_name_full(chain_id, res_no,
                                                        ins_code, "CA",
                                                        alt_conf)
#            rotate_y_scene(10, 0.3) # n_frames frame_interval(degrees)
            with NoBackups(imol):
                with AutoAccept():
                    sphere_refine(radius=refine_residue_sphere_radius, expand=False)
#            rotate_y_scene(10, 0.3)

    set_refine_ramachandran_angles(current_rama_state)
    set_dragged_refinement_steps_per_frame(current_steps_per_frame)


def fix_backbone_range(rng,rng0=0,imol=0,chain='A'):
  resrange=range(rng0,rng)
  for res in resrange:
    for atm in ('CA','N','O','C','H','HA'):
      atm_formated=' {:3}'.format(atm)
      imol=0
      s=[chain,res,'',atm_formated,'']
      #print(s)
      mark_atom_as_fixed(imol,s,1) #1 means fix, 0 means unfix

def unfix_backbone_range(rng,rng0=0,imol=0,chain='A'):
  resrange=range(rng0,rng)
  for res in resrange:
    for atm in ('CA','N','O','C','H','HA'):
      atm_formated=' {:3}'.format(atm)
      imol=0
      s=[chain,res,'',atm_formated,'']
      #print(s)
      mark_atom_as_fixed(imol,s,0) #1 means fix, 0 means unfix

def fix_bb(fix=1):
  imol=active_residue()[0]
  res=active_residue()[1:]
  for atm in ('CA','N','O','C','H','HA'):
    atm_formated=' {:3}'.format(atm)
    res[3]=atm_formated
    mark_atom_as_fixed(imol,res,fix) #1 means fix, 0 means unfix

def fix_res(fix=1):
  imol=active_residue()[0]
  res=active_residue()[2]
  ch=active_residue()[1]

  atoms= residue_info(imol, ch, res, "")
  for a in atoms:
          s=[ch,res,'',a[0][0],'']
          mark_atom_as_fixed(imol,s,1) #1 means fix, 0 means unfix


def unfix_bb():
  fix_bb(fix=0)


def nearby_res(imol=0,rad=2):
  return residues_near_residue(imol,active_residue()[1:4],rad)

def active_res_info():
    active_atom = active_residue()
    imol      = active_atom[0]
    chain_id  = active_atom[1]
    res_no    = active_atom[2]
    ins_code  = active_atom[3]
    atom_name = active_atom[4]
    alt_conf  = active_atom[5]
    print(active_atom)

def refine_with_Nglycan(use_map=True,res_range=3):
    from types import ListType
    active_atom = active_residue()
    if (not active_atom):
        add_status_bar_text("No active residue")
    else:
        if (use_map and not valid_map_molecule_qm(imol_refinement_map())):
            show_select_map_dialog()
        else:
            imol      = active_atom[0]
            chain_id  = active_atom[1]
            res_no    = active_atom[2]
            ins_code  = active_atom[3]
            atom_name = active_atom[4]
            alt_conf  = active_atom[5]
            centred_residue = active_atom[1:4]

            s = []
            aas=range(int(res_no)-res_range,int(res_no)+res_range+1)

            for res in aas: #need to test ends

                s.append([chain_id,res,''])

            potential_nags=residues_near_residue(imol, centred_residue, 2.5)

            for oo in potential_nags:
                if residue_has_hetatms(oo[0],oo[1],oo[2],oo[3])==1:
                    s.append(oo[1:4])
            if s:
              if use_map:
                refine_residues(imol, s)
              else:
                regularize_residues(imol, all_residues)
def sel(chain,res):
    return [[chain,i,''] for i in res]

def sel_range(chain,start,stop):
    return [[chain,i,''] for i in range(start,stop+1)]

def refine_linked(use_map=True,fixbb=True):
    from types import ListType
    active_atom = active_residue()
    if (not active_atom):
        add_status_bar_text("No active residue")
    else:
        if (use_map and not valid_map_molecule_qm(imol_refinement_map())):
            show_select_map_dialog()
        else:
            imol      = active_atom[0]
            chain_id  = active_atom[1]
            res_no    = active_atom[2]
            ins_code  = active_atom[3]
            atom_name = active_atom[4]
            alt_conf  = active_atom[5]
            centred_residue = active_atom[1:4]

            other_residues = residues_near_position(imol,rotation_center(),2.5)
            #residues_near_residue(imol, centred_residue, 2.5)

            other=None
            #print(centred_residue,"added")

            for oo in other_residues:
              o=oo[1:4]
              if o[0] == centred_residue[0]: #same chain
              #  if o[1] != centred_residue[1]+1:
              #    if o[1] != centred_residue[1]-1:
                    other=o
                    print(other,"added")
              else: print(o,"rejected")
            if other is not None:
              all_residues = [centred_residue,other]
              if fixbb:
                for o in all_residues:
                  print("fixing ",o[0],o[1])
                  fix_backbone_range((o[1]+1),rng0=o[1],imol=imol,chain=o[0])


              print "imol: %s residues: %s" %(imol, all_residues)
              if use_map:
                refine_residues(imol, all_residues, "",
                                                         '#soft-mode/hard-mode', False, False)
              else:
                regularize_residues(imol, all_residues)

# Sphere refinement (around radius)


def fix_arg(imol=0,ch='A',rng0=1):
  for res in range(rng0,1034):

    resname=residue_name(imol, ch, res, "")
    '''[
     [[' N  ', ''], [1.0, 68.58, ' N', ''], [190.138, 157.113, 153.696], 127],
     [[' CA ', ''], [1.0, 68.58, ' C', ''], [189.61, 156.65, 154.978], 128], ...]'''
    if resname== 'ARG':

      atoms= residue_info(imol, ch, res, "")
      for a in atoms:
        if a[0][0] != ' NH1':
          s=[ch,res,'',a[0][0],'']
          mark_atom_as_fixed(imol,s,1) #1 means fix, 0 means unfix

      regularize_residues(imol,[[ch,res,'']])
      accept_regularizement()
