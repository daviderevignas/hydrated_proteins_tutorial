import numpy as np
import os
from numba import njit


    
def load_dictionary_loc_density(npz_loaded_file):
    
    box_size=int(npz_loaded_file['box_size'])
    nativ_neigh_list=npz_loaded_file['nativ_neigh_list']
    all_x_pos=npz_loaded_file['all_x_pos']
    all_y_pos=npz_loaded_file['all_y_pos']
    Rg_to_use=npz_loaded_file['Rg_to_use']
    aa_present=npz_loaded_file['aa_present']
    Nc_present=npz_loaded_file['Nc_present']
    Ic_present=npz_loaded_file['Ic_present']
    hw_present=npz_loaded_file['hw_present']
    bw_present=npz_loaded_file['bw_present']
    bulk_HB_present=npz_loaded_file['bulk_HB_present']
    pho_HB_present=npz_loaded_file['pho_HB_present']
    phi_HB_present=npz_loaded_file['phi_HB_present']
    mix_HB_present=npz_loaded_file['mix_HB_present']
    coop_bonds_present=npz_loaded_file['coop_bonds_present']
    coop_bonds_PHO_present=npz_loaded_file['coop_bonds_PHO_present']
    aa_counts=npz_loaded_file['aa_counts']
    Nc_counts=npz_loaded_file['Nc_counts']
    Ic_counts=npz_loaded_file['Ic_counts']
    hw_counts=npz_loaded_file['hw_counts']
    bw_counts=npz_loaded_file['bw_counts']
    hb_b_counts=npz_loaded_file['hb_b_counts']
    hb_pho_counts=npz_loaded_file['hb_pho_counts']
    hb_phi_counts=npz_loaded_file['hb_phi_counts']
    hb_mix_counts=npz_loaded_file['hb_mix_counts']
    coop_b_counts=npz_loaded_file['coop_b_counts']
    coop_pho_counts=npz_loaded_file['coop_pho_counts']
    prot_counts=npz_loaded_file['prot_counts']
    n_cells_tested=npz_loaded_file['n_cells_tested']
    raw_data_folding=npz_loaded_file['raw_data_folding']
    skip_times=int(npz_loaded_file['skip_times'])
    list_of_times_skipped=npz_loaded_file['list_of_times_skipped']
    
    liq_cluster_label=npz_loaded_file['liq_cluster_label']
    sol_cluster_label=npz_loaded_file['sol_cluster_label']
    droplet_cluster_label=npz_loaded_file['droplet_cluster_label']

    prot_belongs_to_droplet_cluster=npz_loaded_file['prot_belongs_to_droplet_cluster']
    droplet_cluster_size_prot=npz_loaded_file['droplet_cluster_size_prot']
    aa_belongs_to_droplet_cluster=npz_loaded_file['aa_belongs_to_droplet_cluster']
    droplet_cluster_size_aa=npz_loaded_file['droplet_cluster_size_aa']


    droplet_cluster_hor_max_size=npz_loaded_file['droplet_cluster_hor_max_size']
    droplet_cluster_ver_max_size=npz_loaded_file['droplet_cluster_ver_max_size']
    all_Nc=npz_loaded_file['all_Nc']
    all_Ic=npz_loaded_file['all_Ic']
    
    cluster_hor_max_size=npz_loaded_file['cluster_hor_max_size']
    cluster_ver_max_size=npz_loaded_file['cluster_ver_max_size']
    sol_cluster_hor_max_size=npz_loaded_file['sol_cluster_hor_max_size']
    sol_cluster_ver_max_size=npz_loaded_file['sol_cluster_ver_max_size']
    
        
    return box_size,nativ_neigh_list,\
        all_x_pos,all_y_pos,Rg_to_use,aa_present,\
        Nc_present,Ic_present,hw_present,bw_present,\
        bulk_HB_present,pho_HB_present,phi_HB_present,\
        mix_HB_present,coop_bonds_present,coop_bonds_PHO_present,\
        aa_counts,Nc_counts,Ic_counts,hw_counts,bw_counts,hb_b_counts,hb_pho_counts,\
        hb_phi_counts,hb_mix_counts,coop_b_counts,coop_pho_counts,prot_counts,n_cells_tested,raw_data_folding,\
        skip_times,list_of_times_skipped,\
        liq_cluster_label,sol_cluster_label,\
        droplet_cluster_label,droplet_cluster_hor_max_size,droplet_cluster_ver_max_size,\
        prot_belongs_to_droplet_cluster,droplet_cluster_size_prot,aa_belongs_to_droplet_cluster,droplet_cluster_size_aa,\
        all_Nc,all_Ic,cluster_hor_max_size,cluster_ver_max_size,sol_cluster_hor_max_size,sol_cluster_ver_max_size
    

        


def build_result_dict_loc_density(folder_with_dicts,selected_Np,selected_T,selected_bs,selected_st):
    
    # CAREFUL! here we make assumptions on the thermodynamic conditions!
    Press=0.0007
    compress_phob=4
    vHB_bulk=.5
    vHB_phob=2-compress_phob*Press
    vHB_phil=.5
    vHB_mix=0.5*vHB_phob+0.5*vHB_phil
   
    
    all_res_dict={}
    if len(selected_Np)==0: selected_Np=[ii for ii in range(1000)]
    if len(selected_T)==0: selected_T=[ii/1000 for ii in range(10000)]
    if len(selected_bs)==0: selected_bs=[ii for ii in range(200)]
    if len(selected_st)==0: selected_st=[ii for ii in range(5)]
    for f in os.listdir(folder_with_dicts):
        if '.npz' in f:
            now_Np=int(f.split('Np_')[-1].split('_')[0])
            now_Tunit=float(f.split('T_')[-1].split('_')[0])
            now_Tdec=float(f.split('T_')[-1].split('_')[1])
            now_T=now_Tunit+now_Tdec/1000
            now_bs=int(f.split('bs_')[-1].split('_')[0])
            now_st=int(f.split('sc_')[-1].split('.')[0])
            if now_Np in selected_Np and\
                now_T in selected_T and\
                now_bs in selected_bs and\
                now_st in selected_st:
                    
                box_size,nativ_neigh_list,\
                    all_x_pos,all_y_pos,Rg_to_use,aa_present,\
                    Nc_present,Ic_present,hw_present,bw_present,\
                    bulk_HB_present,pho_HB_present,phi_HB_present,\
                    mix_HB_present,coop_bonds_present,coop_bonds_PHO_present,\
                    aa_counts,Nc_counts,Ic_counts,hw_counts,bw_counts,hb_b_counts,hb_pho_counts,\
                    hb_phi_counts,hb_mix_counts,coop_b_counts,coop_pho_counts,prot_counts,n_cells_tested,raw_data_folding,\
                    skip_times,list_of_times_skipped,\
                    liq_cluster_label,sol_cluster_label,\
                    droplet_cluster_label,droplet_cluster_hor_max_size,droplet_cluster_ver_max_size,\
                    prot_belongs_to_droplet_cluster,droplet_cluster_size_prot,aa_belongs_to_droplet_cluster,droplet_cluster_size_aa,\
                    all_Nc,all_Ic,\
                    cluster_hor_max_size,cluster_ver_max_size,\
                    sol_cluster_hor_max_size,sol_cluster_ver_max_size,\
                    =load_dictionary_loc_density(np.load(folder_with_dicts+f))
                
                # liq_cluster_size,prot_belongs_to_liq_cluster,\
                #     sol_cluster_size,prot_belongs_to_sol_cluster,\
                #     tot_liq_cluster_sizes,tot_liq_cluster_volumes,tot_liq_cluster_volume_fractions,\
                #     tot_sol_cluster_sizes,tot_sol_cluster_volumes,tot_sol_cluster_volume_fractions,\
                #     cluster_hor_max_size,cluster_ver_max_size,\
                #     sol_cluster_hor_max_size,sol_cluster_ver_max_size,\
                #     liq_cluster_gyr_radius,sol_cluster_gyr_radius
                V_cell_proper=(
                        raw_data_folding[:,3,:]
                        -raw_data_folding[:,4,:]*vHB_bulk
                        -raw_data_folding[:,5,:]*vHB_phob
                        -raw_data_folding[:,6,:]*vHB_phil
                        -raw_data_folding[:,7,:]*vHB_mix
                    )/(box_size**2)
                V_probe_hb=.5*(hb_b_counts*vHB_bulk+hb_pho_counts*vHB_phob+hb_phi_counts*vHB_phil+hb_mix_counts*vHB_mix)
                V_probe=V_probe_hb+V_cell_proper*n_cells_tested
                V_cell=V_cell_proper+.5*(bulk_HB_present*vHB_bulk+pho_HB_present*vHB_phob+phi_HB_present*vHB_phil+mix_HB_present*vHB_mix)
                
                all_res_dict[(
                        now_Np,now_T,now_bs,now_st
                        )]= \
                    box_size,nativ_neigh_list,\
                    all_x_pos,all_y_pos,Rg_to_use,aa_present,\
                    Nc_present,Ic_present,hw_present,bw_present,\
                    bulk_HB_present,pho_HB_present,phi_HB_present,\
                    mix_HB_present,coop_bonds_present,coop_bonds_PHO_present,\
                    aa_counts,Nc_counts,Ic_counts,hw_counts,bw_counts,hb_b_counts,hb_pho_counts,\
                    hb_phi_counts,hb_mix_counts,coop_b_counts,coop_pho_counts,prot_counts,n_cells_tested,raw_data_folding,\
                    V_cell_proper,V_probe_hb,V_probe,V_cell,\
                    skip_times,list_of_times_skipped,\
                    liq_cluster_label,sol_cluster_label,\
                    droplet_cluster_label,droplet_cluster_hor_max_size,droplet_cluster_ver_max_size,\
                    cluster_hor_max_size,cluster_ver_max_size,\
                    sol_cluster_hor_max_size,sol_cluster_ver_max_size,\
                    all_Nc,all_Ic,\
                    cluster_hor_max_size,cluster_ver_max_size,\
                    sol_cluster_hor_max_size,sol_cluster_ver_max_size,\
                    prot_belongs_to_droplet_cluster,droplet_cluster_size_prot,aa_belongs_to_droplet_cluster,droplet_cluster_size_aa
                    # liq_cluster_size,prot_belongs_to_liq_cluster,\
                    # sol_cluster_size,prot_belongs_to_sol_cluster,\
                    # tot_liq_cluster_sizes,tot_liq_cluster_volumes,tot_liq_cluster_volume_fractions,\
                    # tot_sol_cluster_sizes,tot_sol_cluster_volumes,tot_sol_cluster_volume_fractions,\
                    # cluster_hor_max_size,cluster_ver_max_size,\
                    # sol_cluster_hor_max_size,sol_cluster_ver_max_size,\
                    # liq_cluster_gyr_radius,sol_cluster_gyr_radius,\
                print(f,now_Np,now_T,now_bs,now_st)
    return all_res_dict






# %%
def build_result_dict(folder_with_dicts,selected_Np,selected_T,selected_bs,selected_st):
    all_res_dict={}
    if len(selected_Np)==0: selected_Np=[ii for ii in range(1000)]
    if len(selected_T)==0: selected_T=[ii/1000 for ii in range(10000)]
    if len(selected_bs)==0: selected_bs=[ii for ii in range(200)]
    if len(selected_st)==0: selected_st=[ii for ii in range(5)]
    for f in os.listdir(folder_with_dicts):
        if '.npz' in f:
            now_Np=int(f.split('Np_')[-1].split('_')[0])
            now_Tunit=float(f.split('T_')[-1].split('_')[0])
            now_Tdec=float(f.split('T_')[-1].split('_')[1])
            now_T=now_Tunit+now_Tdec/1000
            now_bs=int(f.split('bs_')[-1].split('_')[0])
            now_st=int(f.split('sc_')[-1].split('.')[0])
            if now_Np in selected_Np and\
                now_T in selected_T and\
                now_bs in selected_bs and\
                now_st in selected_st:
                    
                box_size,nativ_neigh_list,path_to_array_replica_folder,\
                    all_Nc,all_Ic,list_of_times,\
                    raw_data_folding,\
                    n_monomers,n_proteins,n_snapshots,n_replicas,n_native_bonds_per_protein,\
                    liq_cluster_size,prot_belongs_to_liq_cluster,\
                    sol_cluster_size,prot_belongs_to_sol_cluster,\
                    tot_liq_cluster_sizes,tot_liq_cluster_volumes,tot_liq_cluster_volume_fractions,\
                    tot_sol_cluster_sizes,tot_sol_cluster_volumes,tot_sol_cluster_volume_fractions,\
                    cluster_hor_max_size,cluster_ver_max_size,\
                    sol_cluster_hor_max_size,sol_cluster_ver_max_size,\
                    liq_cluster_gyr_radius,sol_cluster_gyr_radius,\
                    MC_V_probe,probed_cells,\
                    MC_V_a,MC_V_h,MC_V_b,\
                    MC_N_a,MC_N_h,MC_N_b,\
                    MC_diff_prot_counts,MC_prot_wat_energy,MC_prot_energy,\
                    MC_N_hb_b,MC_N_hb_phi,MC_N_hb_pho,MC_N_hb_mix,\
                    MC_N_coop_b,MC_N_coop_pho,MC_N_coop_phi,MC_N_coop_mix,\
                    MC_hb_b_energy,MC_hb_phi_energy,MC_hb_pho_energy,MC_hb_mix_energy,\
                    MC_coop_b_energy,MC_coop_pho_energy,MC_coop_phi_energy,MC_coop_mix_energy,\
                    MC_LJ_energy,\
                    MC_N_sig1_h,MC_N_sig2_h,MC_N_sig3_h,MC_N_sig4_h,MC_N_sig5_h,MC_N_sig6_h,\
                    MC_N_sig1_b,MC_N_sig2_b,MC_N_sig3_b,MC_N_sig4_b,MC_N_sig5_b,MC_N_sig6_b,\
                    MC_N_hb_b_max_possible,MC_N_hb_h_max_possible=load_dictionary_all_simple(np.load(folder_with_dicts+f))
                all_res_dict[(
                        now_Np,now_T,now_bs,now_st
                        )]= \
                        box_size,nativ_neigh_list,path_to_array_replica_folder,\
                        all_Nc,all_Ic,list_of_times,\
                        raw_data_folding,\
                        n_monomers,n_proteins,n_snapshots,n_replicas,n_native_bonds_per_protein,\
                        liq_cluster_size,prot_belongs_to_liq_cluster,\
                        sol_cluster_size,prot_belongs_to_sol_cluster,\
                        tot_liq_cluster_sizes,tot_liq_cluster_volumes,tot_liq_cluster_volume_fractions,\
                        tot_sol_cluster_sizes,tot_sol_cluster_volumes,tot_sol_cluster_volume_fractions,\
                        cluster_hor_max_size,cluster_ver_max_size,\
                        sol_cluster_hor_max_size,sol_cluster_ver_max_size,\
                        liq_cluster_gyr_radius,sol_cluster_gyr_radius,\
                        MC_V_probe,probed_cells,\
                        MC_V_a,MC_V_h,MC_V_b,\
                        MC_N_a,MC_N_h,MC_N_b,\
                        MC_diff_prot_counts,MC_prot_wat_energy,MC_prot_energy,\
                        MC_N_hb_b,MC_N_hb_phi,MC_N_hb_pho,MC_N_hb_mix,\
                        MC_N_coop_b,MC_N_coop_pho,MC_N_coop_phi,MC_N_coop_mix,\
                        MC_hb_b_energy,MC_hb_phi_energy,MC_hb_pho_energy,MC_hb_mix_energy,\
                        MC_coop_b_energy,MC_coop_pho_energy,MC_coop_phi_energy,MC_coop_mix_energy,\
                        MC_LJ_energy,\
                        MC_N_sig1_h,MC_N_sig2_h,MC_N_sig3_h,MC_N_sig4_h,MC_N_sig5_h,MC_N_sig6_h,\
                        MC_N_sig1_b,MC_N_sig2_b,MC_N_sig3_b,MC_N_sig4_b,MC_N_sig5_b,MC_N_sig6_b,\
                        MC_N_hb_b_max_possible,MC_N_hb_h_max_possible
                
                print(f,now_Np,now_T,now_bs,now_st)
    return all_res_dict

# %%

def load_dictionary_all_simple(npz_loaded_file):
        
    box_size=npz_loaded_file['box_size']
    nativ_neigh_list=npz_loaded_file['nativ_neigh_list']
    path_to_array_replica_folder=npz_loaded_file['path_to_array_replica_folder']
    all_Nc=npz_loaded_file['all_Nc']
    all_Ic=npz_loaded_file['all_Ic']
    list_of_times=npz_loaded_file['list_of_times']
    raw_data_folding=npz_loaded_file['raw_data_folding']
    n_monomers=npz_loaded_file['n_monomers']
    n_proteins=npz_loaded_file['n_proteins']
    n_snapshots=npz_loaded_file['n_snapshots']
    n_replicas=npz_loaded_file['n_replicas']
    n_native_bonds_per_protein=npz_loaded_file['n_native_bonds_per_protein']
    liq_cluster_size=npz_loaded_file['liq_cluster_size']
    prot_belongs_to_liq_cluster=npz_loaded_file['prot_belongs_to_liq_cluster']
    sol_cluster_size=npz_loaded_file['sol_cluster_size']
    prot_belongs_to_sol_cluster=npz_loaded_file['prot_belongs_to_sol_cluster']
    tot_liq_cluster_sizes=npz_loaded_file['tot_liq_cluster_sizes']
    tot_liq_cluster_volumes=npz_loaded_file['tot_liq_cluster_volumes']
    tot_liq_cluster_volume_fractions=npz_loaded_file['tot_liq_cluster_volume_fractions']
    tot_sol_cluster_sizes=npz_loaded_file['tot_sol_cluster_sizes']
    tot_sol_cluster_volumes=npz_loaded_file['tot_sol_cluster_volumes']
    tot_sol_cluster_volume_fractions=npz_loaded_file['tot_sol_cluster_volume_fractions']
    cluster_hor_max_size=npz_loaded_file['cluster_hor_max_size']
    cluster_ver_max_size=npz_loaded_file['cluster_ver_max_size']
    sol_cluster_hor_max_size=npz_loaded_file['sol_cluster_hor_max_size']
    sol_cluster_ver_max_size=npz_loaded_file['sol_cluster_ver_max_size']
    liq_cluster_gyr_radius=npz_loaded_file['liq_cluster_gyr_radius']
    sol_cluster_gyr_radius=npz_loaded_file['sol_cluster_gyr_radius']
    MC_V_probe=npz_loaded_file['MC_V_probe']
    probed_cells=npz_loaded_file['probed_cells']
    MC_V_a=npz_loaded_file['MC_V_a']
    MC_V_h=npz_loaded_file['MC_V_h']
    MC_V_b=npz_loaded_file['MC_V_b']
    MC_N_a=npz_loaded_file['MC_N_a']
    MC_N_h=npz_loaded_file['MC_N_h']
    MC_N_b=npz_loaded_file['MC_N_b']
    MC_diff_prot_counts=npz_loaded_file['MC_diff_prot_counts']
    MC_prot_wat_energy=npz_loaded_file['MC_prot_wat_energy']
    MC_prot_energy=npz_loaded_file['MC_prot_energy']
    MC_N_hb_b=npz_loaded_file['MC_N_hb_b']
    MC_N_hb_phi=npz_loaded_file['MC_N_hb_phi']
    MC_N_hb_pho=npz_loaded_file['MC_N_hb_pho']
    MC_N_hb_mix=npz_loaded_file['MC_N_hb_mix']
    MC_N_coop_b=npz_loaded_file['MC_N_coop_b']
    MC_N_coop_pho=npz_loaded_file['MC_N_coop_pho']
    MC_N_coop_phi=npz_loaded_file['MC_N_coop_phi']
    MC_N_coop_mix=npz_loaded_file['MC_N_coop_mix']
    MC_hb_b_energy=npz_loaded_file['MC_hb_b_energy']
    MC_hb_phi_energy=npz_loaded_file['MC_hb_phi_energy']
    MC_hb_pho_energy=npz_loaded_file['MC_hb_pho_energy']
    MC_hb_mix_energy=npz_loaded_file['MC_hb_mix_energy']
    MC_coop_b_energy=npz_loaded_file['MC_coop_b_energy']
    MC_coop_pho_energy=npz_loaded_file['MC_coop_pho_energy']
    MC_coop_phi_energy=npz_loaded_file['MC_coop_phi_energy']
    MC_coop_mix_energy=npz_loaded_file['MC_coop_mix_energy']
    MC_LJ_energy=npz_loaded_file['MC_LJ_energy']
    MC_N_sig1_h=npz_loaded_file['MC_N_sig1_h']
    MC_N_sig2_h=npz_loaded_file['MC_N_sig2_h']
    MC_N_sig3_h=npz_loaded_file['MC_N_sig3_h']
    MC_N_sig4_h=npz_loaded_file['MC_N_sig4_h']
    MC_N_sig5_h=npz_loaded_file['MC_N_sig5_h']
    MC_N_sig6_h=npz_loaded_file['MC_N_sig6_h']
    MC_N_sig1_b=npz_loaded_file['MC_N_sig1_b']
    MC_N_sig2_b=npz_loaded_file['MC_N_sig2_b']
    MC_N_sig3_b=npz_loaded_file['MC_N_sig3_b']
    MC_N_sig4_b=npz_loaded_file['MC_N_sig4_b']
    MC_N_sig5_b=npz_loaded_file['MC_N_sig5_b']
    MC_N_sig6_b=npz_loaded_file['MC_N_sig6_b']
    MC_N_hb_b_max_possible=npz_loaded_file['MC_N_hb_b_max_possible']
    MC_N_hb_h_max_possible=npz_loaded_file['MC_N_hb_h_max_possible']
    
    return box_size,nativ_neigh_list,path_to_array_replica_folder,\
        all_Nc,all_Ic,list_of_times,\
        raw_data_folding,\
        n_monomers,n_proteins,n_snapshots,n_replicas,n_native_bonds_per_protein,\
        liq_cluster_size,prot_belongs_to_liq_cluster,\
        sol_cluster_size,prot_belongs_to_sol_cluster,\
        tot_liq_cluster_sizes,tot_liq_cluster_volumes,tot_liq_cluster_volume_fractions,\
        tot_sol_cluster_sizes,tot_sol_cluster_volumes,tot_sol_cluster_volume_fractions,\
        cluster_hor_max_size,cluster_ver_max_size,\
        sol_cluster_hor_max_size,sol_cluster_ver_max_size,\
        liq_cluster_gyr_radius,sol_cluster_gyr_radius,\
        MC_V_probe,probed_cells,\
        MC_V_a,MC_V_h,MC_V_b,\
        MC_N_a,MC_N_h,MC_N_b,\
        MC_diff_prot_counts,MC_prot_wat_energy,MC_prot_energy,\
        MC_N_hb_b,MC_N_hb_phi,MC_N_hb_pho,MC_N_hb_mix,\
        MC_N_coop_b,MC_N_coop_pho,MC_N_coop_phi,MC_N_coop_mix,\
        MC_hb_b_energy,MC_hb_phi_energy,MC_hb_pho_energy,MC_hb_mix_energy,\
        MC_coop_b_energy,MC_coop_pho_energy,MC_coop_phi_energy,MC_coop_mix_energy,\
        MC_LJ_energy,\
        MC_N_sig1_h,MC_N_sig2_h,MC_N_sig3_h,MC_N_sig4_h,MC_N_sig5_h,MC_N_sig6_h,\
        MC_N_sig1_b,MC_N_sig2_b,MC_N_sig3_b,MC_N_sig4_b,MC_N_sig5_b,MC_N_sig6_b,\
        MC_N_hb_b_max_possible,MC_N_hb_h_max_possible






# %%
def unwrap_MC_results(npz_loaded_file_list):
    if type(npz_loaded_file_list) == list:
        MC_V_probe=npz_loaded_file_list[0]['MC_V_probe']
        probed_cells=int(npz_loaded_file_list[0]['probed_cells'])
        MC_V_a=npz_loaded_file_list[0]['MC_V_a']
        MC_V_h=npz_loaded_file_list[0]['MC_V_h']
        MC_V_b=npz_loaded_file_list[0]['MC_V_b']
        MC_N_a=npz_loaded_file_list[0]['MC_N_a']
        MC_N_h=npz_loaded_file_list[0]['MC_N_h']
        MC_N_b=npz_loaded_file_list[0]['MC_N_b']
        MC_diff_prot_counts=npz_loaded_file_list[0]['MC_diff_prot_counts']
        MC_prot_wat_energy=npz_loaded_file_list[0]['MC_prot_wat_energy']
        MC_prot_energy=npz_loaded_file_list[0]['MC_prot_energy']
        MC_N_hb_b=npz_loaded_file_list[0]['MC_N_hb_b']
        MC_N_hb_phi=npz_loaded_file_list[0]['MC_N_hb_phi']
        MC_N_hb_pho=npz_loaded_file_list[0]['MC_N_hb_pho']
        MC_N_hb_mix=npz_loaded_file_list[0]['MC_N_hb_mix']
        MC_N_coop_b=npz_loaded_file_list[0]['MC_N_coop_b']
        MC_N_coop_pho=npz_loaded_file_list[0]['MC_N_coop_pho']
        MC_N_coop_phi=npz_loaded_file_list[0]['MC_N_coop_phi']
        MC_N_coop_mix=npz_loaded_file_list[0]['MC_N_coop_mix']
        MC_hb_b_energy=npz_loaded_file_list[0]['MC_hb_b_energy']
        MC_hb_phi_energy=npz_loaded_file_list[0]['MC_hb_phi_energy']
        MC_hb_pho_energy=npz_loaded_file_list[0]['MC_hb_pho_energy']
        MC_hb_mix_energy=npz_loaded_file_list[0]['MC_hb_mix_energy']
        MC_coop_b_energy=npz_loaded_file_list[0]['MC_coop_b_energy']
        MC_coop_pho_energy=npz_loaded_file_list[0]['MC_coop_pho_energy']
        MC_coop_phi_energy=npz_loaded_file_list[0]['MC_coop_phi_energy']
        MC_coop_mix_energy=npz_loaded_file_list[0]['MC_coop_mix_energy']
        MC_LJ_energy=npz_loaded_file_list[0]['MC_LJ_energy']
        MC_N_sig1_h=npz_loaded_file_list[0]['MC_N_sig1_h']
        MC_N_sig2_h=npz_loaded_file_list[0]['MC_N_sig2_h']
        MC_N_sig3_h=npz_loaded_file_list[0]['MC_N_sig3_h']
        MC_N_sig4_h=npz_loaded_file_list[0]['MC_N_sig4_h']
        MC_N_sig5_h=npz_loaded_file_list[0]['MC_N_sig5_h']
        MC_N_sig6_h=npz_loaded_file_list[0]['MC_N_sig6_h']
        MC_N_sig1_b=npz_loaded_file_list[0]['MC_N_sig1_b']
        MC_N_sig2_b=npz_loaded_file_list[0]['MC_N_sig2_b']
        MC_N_sig3_b=npz_loaded_file_list[0]['MC_N_sig3_b']
        MC_N_sig4_b=npz_loaded_file_list[0]['MC_N_sig4_b']
        MC_N_sig5_b=npz_loaded_file_list[0]['MC_N_sig5_b']
        MC_N_sig6_b=npz_loaded_file_list[0]['MC_N_sig6_b']
        MC_N_hb_b_max_possible=npz_loaded_file_list[0]['MC_N_hb_b_max_possible']
        MC_N_hb_h_max_possible=npz_loaded_file_list[0]['MC_N_hb_h_max_possible']
    for i_file in range(1,len(npz_loaded_file_list)):
        MC_V_probe=np.append(MC_V_probe,npz_loaded_file_list[i_file]['MC_V_probe'],axis=0)
        MC_V_a=np.append(MC_V_a,npz_loaded_file_list[i_file]['MC_V_a'],axis=0)
        MC_V_h=np.append(MC_V_h,npz_loaded_file_list[i_file]['MC_V_h'],axis=0)
        MC_V_b=np.append(MC_V_b,npz_loaded_file_list[i_file]['MC_V_b'],axis=0)
        MC_N_a=np.append(MC_N_a,npz_loaded_file_list[i_file]['MC_N_a'],axis=0)
        MC_N_h=np.append(MC_N_h,npz_loaded_file_list[i_file]['MC_N_h'],axis=0)
        MC_N_b=np.append(MC_N_b,npz_loaded_file_list[i_file]['MC_N_b'],axis=0)
        MC_diff_prot_counts=np.append(MC_diff_prot_counts,npz_loaded_file_list[i_file]['MC_diff_prot_counts'],axis=0)
        MC_prot_wat_energy=np.append(MC_prot_wat_energy,npz_loaded_file_list[i_file]['MC_prot_wat_energy'],axis=0)
        MC_prot_energy=np.append(MC_prot_energy,npz_loaded_file_list[i_file]['MC_prot_energy'],axis=0)
        MC_N_hb_b=np.append(MC_N_hb_b,npz_loaded_file_list[i_file]['MC_N_hb_b'],axis=0)
        MC_N_hb_phi=np.append(MC_N_hb_phi,npz_loaded_file_list[i_file]['MC_N_hb_phi'],axis=0)
        MC_N_hb_pho=np.append(MC_N_hb_pho,npz_loaded_file_list[i_file]['MC_N_hb_pho'],axis=0)
        MC_N_hb_mix=np.append(MC_N_hb_mix,npz_loaded_file_list[i_file]['MC_N_hb_mix'],axis=0)
        MC_N_coop_b=np.append(MC_N_coop_b,npz_loaded_file_list[i_file]['MC_N_coop_b'],axis=0)
        MC_N_coop_pho=np.append(MC_N_coop_pho,npz_loaded_file_list[i_file]['MC_N_coop_pho'],axis=0)
        MC_N_coop_phi=np.append(MC_N_coop_phi,npz_loaded_file_list[i_file]['MC_N_coop_phi'],axis=0)
        MC_N_coop_mix=np.append(MC_N_coop_mix,npz_loaded_file_list[i_file]['MC_N_coop_mix'],axis=0)
        MC_hb_b_energy=np.append(MC_hb_b_energy,npz_loaded_file_list[i_file]['MC_hb_b_energy'],axis=0)
        MC_hb_phi_energy=np.append(MC_hb_phi_energy,npz_loaded_file_list[i_file]['MC_hb_phi_energy'],axis=0)
        MC_hb_pho_energy=np.append(MC_hb_pho_energy,npz_loaded_file_list[i_file]['MC_hb_pho_energy'],axis=0)
        MC_hb_mix_energy=np.append(MC_hb_mix_energy,npz_loaded_file_list[i_file]['MC_hb_mix_energy'],axis=0)
        MC_coop_b_energy=np.append(MC_coop_b_energy,npz_loaded_file_list[i_file]['MC_coop_b_energy'],axis=0)
        MC_coop_pho_energy=np.append(MC_coop_pho_energy,npz_loaded_file_list[i_file]['MC_coop_pho_energy'],axis=0)
        MC_coop_phi_energy=np.append(MC_coop_phi_energy,npz_loaded_file_list[i_file]['MC_coop_phi_energy'],axis=0)
        MC_coop_mix_energy=np.append(MC_coop_mix_energy,npz_loaded_file_list[i_file]['MC_coop_mix_energy'],axis=0)
        MC_LJ_energy=np.append(MC_LJ_energy,npz_loaded_file_list[i_file]['MC_LJ_energy'],axis=0)
        MC_N_sig1_h=np.append(MC_N_sig1_h,npz_loaded_file_list[i_file]['MC_N_sig1_h'],axis=0)
        MC_N_sig2_h=np.append(MC_N_sig2_h,npz_loaded_file_list[i_file]['MC_N_sig2_h'],axis=0)
        MC_N_sig3_h=np.append(MC_N_sig3_h,npz_loaded_file_list[i_file]['MC_N_sig3_h'],axis=0)
        MC_N_sig4_h=np.append(MC_N_sig4_h,npz_loaded_file_list[i_file]['MC_N_sig4_h'],axis=0)
        MC_N_sig5_h=np.append(MC_N_sig5_h,npz_loaded_file_list[i_file]['MC_N_sig5_h'],axis=0)
        MC_N_sig6_h=np.append(MC_N_sig6_h,npz_loaded_file_list[i_file]['MC_N_sig6_h'],axis=0)
        MC_N_sig1_b=np.append(MC_N_sig1_b,npz_loaded_file_list[i_file]['MC_N_sig1_b'],axis=0)
        MC_N_sig2_b=np.append(MC_N_sig2_b,npz_loaded_file_list[i_file]['MC_N_sig2_b'],axis=0)
        MC_N_sig3_b=np.append(MC_N_sig3_b,npz_loaded_file_list[i_file]['MC_N_sig3_b'],axis=0)
        MC_N_sig4_b=np.append(MC_N_sig4_b,npz_loaded_file_list[i_file]['MC_N_sig4_b'],axis=0)
        MC_N_sig5_b=np.append(MC_N_sig5_b,npz_loaded_file_list[i_file]['MC_N_sig5_b'],axis=0)
        MC_N_sig6_b=np.append(MC_N_sig6_b,npz_loaded_file_list[i_file]['MC_N_sig6_b'],axis=0)
        MC_N_hb_b_max_possible=np.append(MC_N_hb_b_max_possible,npz_loaded_file_list[i_file]['MC_N_hb_b_max_possible'],axis=0)
        MC_N_hb_h_max_possible=np.append(MC_N_hb_h_max_possible,npz_loaded_file_list[i_file]['MC_N_hb_h_max_possible'],axis=0)
    return MC_V_probe,probed_cells,\
                MC_V_a,MC_V_h,MC_V_b,\
                MC_N_a,MC_N_h,MC_N_b,\
                MC_diff_prot_counts,MC_prot_wat_energy,MC_prot_energy,\
                MC_N_hb_b,MC_N_hb_phi,MC_N_hb_pho,MC_N_hb_mix,\
                MC_N_coop_b,MC_N_coop_pho,MC_N_coop_phi,MC_N_coop_mix,\
                MC_hb_b_energy,MC_hb_phi_energy,MC_hb_pho_energy,MC_hb_mix_energy,\
                MC_coop_b_energy,MC_coop_pho_energy,MC_coop_phi_energy,MC_coop_mix_energy,\
                MC_LJ_energy,\
                MC_N_sig1_h,MC_N_sig2_h,MC_N_sig3_h,MC_N_sig4_h,MC_N_sig5_h,MC_N_sig6_h,\
                MC_N_sig1_b,MC_N_sig2_b,MC_N_sig3_b,MC_N_sig4_b,MC_N_sig5_b,MC_N_sig6_b,\
                MC_N_hb_b_max_possible,MC_N_hb_h_max_possible



def save_analysis_at_folder(path_to_array_replica_folder,n_proteins,n_monomers,box_size):
    
    # ATTENTION, this is declined for A0 protein:
    # nativ_neigh_list=read_target_structures('target_structures.dat',n_monomers)
    nativ_neigh_list=np.array([11, 10,  9,  8,  7, -1, 17, 16, 15, 14, 13, -1, 23, 22, 21, 20, 19,
       -1, 29, 28, 27, 26, 25, -1, 35, 34, 33, 32, 31, -1, -1, -1, -1, -1,
       -1, -1], dtype=np.int32)
    
    # read all protein configurations output from a folder containing replica folders
    # output variables are:
    # -> all_x_pos[i_mon,i_prot,i_t,i_replica]
    # -> all_x_pos[i_mon,i_prot,i_t,i_replica]
    # -> list_of_times[i_t]
    max_n_snapshots=int(1e5)
    n_replicas=0
    replica_indexes=np.array([],dtype=np.int32)
    if 'all_x_pos' in locals(): del(all_x_pos)
    for subdir, dirs, files in os.walk(path_to_array_replica_folder):
        for file in files:
            filepath = subdir + os.sep + file
            if "Protein_Configurations" in filepath:
                n_replicas += 1
                now_replica_index=int(filepath.split('/')[-2])
                replica_indexes= np.append(replica_indexes,now_replica_index)
                raw_conf_data=read_conf_file(filepath,n_proteins,n_monomers,max_n_snapshots)
                x_pos,y_pos,list_of_times = convert_raw_conf_data(raw_conf_data,n_proteins,n_monomers)
                if 'all_x_pos' not in locals():
                    all_x_pos = np.copy(x_pos)
                    all_y_pos = np.copy(y_pos)
                    all_x_pos = np.expand_dims(all_x_pos,axis=3)
                    all_y_pos = np.expand_dims(all_y_pos,axis=3)
                else:
                    all_x_pos = np.concatenate((all_x_pos, np.expand_dims(x_pos,axis=3)), axis=3)
                    all_y_pos = np.concatenate((all_y_pos, np.expand_dims(y_pos,axis=3)), axis=3)
                    
    n_snapshots = len(list_of_times)

    all_x_pos=all_x_pos[:,:,:,np.argsort(replica_indexes)]
    all_y_pos=all_y_pos[:,:,:,np.argsort(replica_indexes)]
    
    
    # read all thermodynamic output from a folder containing replica foldersbox
    #
    # output variables are:
    # -> raw_data_folding[i_t,:,i_replica]

    raw_data_folding=np.zeros((n_snapshots,12,n_replicas),dtype='double')
    i_replica=0
    replica_indexes=np.array([],dtype=np.int32)
    for subdir, dirs, files in os.walk(path_to_array_replica_folder):
        for file in files:
            filepath = subdir + os.sep + file
            if "Data_folding" in filepath:
                now_replica_index=int(filepath.split('/')[-2])
                replica_indexes= np.append(replica_indexes,now_replica_index)
                now_loaded=np.loadtxt(filepath,dtype='double')
                if now_loaded.shape[1]==10:
                    raw_data_folding[:,:,i_replica]=np.concatenate((now_loaded,np.zeros((now_loaded.shape[0],2))),axis=1)
                else:
                    raw_data_folding[:,:,i_replica]=now_loaded
                i_replica=i_replica+1

    raw_data_folding=raw_data_folding[:,:,np.argsort(replica_indexes)]
    
    
    
    all_hor_HB=np.zeros((n_snapshots,box_size**2,n_replicas),dtype=np.int32)
    i_replica=0
    replica_indexes=np.array([],dtype=np.int32)
    for subdir, dirs, files in os.walk(path_to_array_replica_folder):
        for file in files:
            filepath = subdir + os.sep + file
            if "hor_HB" in filepath:
                now_replica_index=int(filepath.split('/')[-2])
                replica_indexes= np.append(replica_indexes,now_replica_index)
                all_hor_HB[:,:,i_replica]=np.loadtxt(filepath,dtype=np.int32)
                i_replica=i_replica+1

    all_hor_HB=all_hor_HB[:,:,np.argsort(replica_indexes)]
    
    all_ver_HB=np.zeros((n_snapshots,box_size**2,n_replicas),dtype=np.int32)
    i_replica=0
    replica_indexes=np.array([],dtype=np.int32)
    for subdir, dirs, files in os.walk(path_to_array_replica_folder):
        for file in files:
            filepath = subdir + os.sep + file
            if "ver_HB" in filepath:
                now_replica_index=int(filepath.split('/')[-2])
                replica_indexes= np.append(replica_indexes,now_replica_index)
                all_ver_HB[:,:,i_replica]=np.loadtxt(filepath,dtype=np.int32)
                i_replica=i_replica+1

    all_ver_HB=all_ver_HB[:,:,np.argsort(replica_indexes)]
    
    all_coop_bonds=np.zeros((n_snapshots,box_size**2,n_replicas),dtype=np.int32)
    i_replica=0
    replica_indexes=np.array([],dtype=np.int32)
    for subdir, dirs, files in os.walk(path_to_array_replica_folder):
        for file in files:
            filepath = subdir + os.sep + file
            if ("coop_bonds" in filepath) and not("PHO" in filepath):
                now_replica_index=int(filepath.split('/')[-2])
                replica_indexes= np.append(replica_indexes,now_replica_index)
                all_coop_bonds[:,:,i_replica]=np.loadtxt(filepath,dtype=np.int32)
                i_replica=i_replica+1

    all_coop_bonds=all_coop_bonds[:,:,np.argsort(replica_indexes)]
    
    all_coop_bonds_PHO=np.zeros((n_snapshots,box_size**2,n_replicas),dtype=np.int32)
    i_replica=0
    replica_indexes=np.array([],dtype=np.int32)
    for subdir, dirs, files in os.walk(path_to_array_replica_folder):
        for file in files:
            filepath = subdir + os.sep + file
            if ("coop_bonds" in filepath) and ("PHO" in filepath):
                now_replica_index=int(filepath.split('/')[-2])
                replica_indexes= np.append(replica_indexes,now_replica_index)
                all_coop_bonds_PHO[:,:,i_replica]=np.loadtxt(filepath,dtype=np.int32)
                i_replica=i_replica+1

    all_coop_bonds_PHO=all_coop_bonds_PHO[:,:,np.argsort(replica_indexes)]
    
    bulk_HB_present,pho_HB_present,phi_HB_present,mix_HB_present,\
    coop_bonds_present,coop_bonds_PHO_present=convert_HB_to_mol_states(all_hor_HB,all_ver_HB,
                                                            all_coop_bonds,all_coop_bonds_PHO,
                                                            box_size)
    

    # data_folding_labels = ["t","system energy","protein energy",
    #                     "system volume","$\#$ H bonds bulk","$\#$ H bonds PHO",
    #                     "$\#$ H bonds PHI","$\#$ H bonds mixed","$\#$ coop bonds",
    #                     "$\#$ coop bonds PHI"]
    

    # analyze all_x_pos and all_y_pos variables to extract 
    # -> native contacts:               all_Nc[i_replica,i_t]
    # -> inter prot contacts:           all_Ic[i_replica,i_t]
    # -> native contacts per prot:      all_Nc_per_prot[i_prot,i_replica,i_t]
    all_Nc=np.zeros((n_replicas,n_snapshots),dtype=np.int32)
    all_Ic=np.zeros((n_replicas,n_snapshots),dtype=np.int32)
    for i_replica in range(n_replicas):
        all_Nc[i_replica,:]=count_all_native_neigh(all_x_pos[:,:,:,i_replica],all_y_pos[:,:,:,i_replica],nativ_neigh_list,box_size)
        all_Ic[i_replica,:]=count_all_inter_contacts(all_x_pos[:,:,:,i_replica],all_y_pos[:,:,:,i_replica],box_size)
    all_Nc_per_prot=count_fold_unfold_prot(all_x_pos,all_y_pos,box_size,nativ_neigh_list)
    
    # analyze all_x_pos and all_y_pos variables to extract    
    aa_present,Nc_present,Ic_present,hw_present,bw_present,\
    aa_counts,Nc_counts,Ic_counts,hw_counts,bw_counts,\
    hb_b_counts,hb_pho_counts,hb_phi_counts,hb_mix_counts,coop_b_counts,coop_pho_counts,\
    bulk_HB_present,pho_HB_present,phi_HB_present,mix_HB_present,coop_bonds_present,coop_bonds_PHO_present,\
    prot_counts,\
    n_cells_tested=calculate_all_local(all_x_pos,all_y_pos,box_size,10,nativ_neigh_list,\
    bulk_HB_present,pho_HB_present,phi_HB_present,mix_HB_present,coop_bonds_present,coop_bonds_PHO_present)
    
    
    
    
    
    st0=np.zeros((n_snapshots,box_size**2,n_replicas),dtype=np.int32)
    n_replicas=0
    replica_indexes=np.array([],dtype=np.int32)
    for subdir, dirs, files in os.walk(path_to_array_replica_folder):
        for file in files:
            filepath = subdir + os.sep + file
            if "status0" in filepath:
                n_replicas += 1
                now_replica_index=int(filepath.split('/')[-2])
                replica_indexes= np.append(replica_indexes,now_replica_index)
                st0[:,:,n_replicas-1]=np.loadtxt(filepath)
    st0=st0[:,:,np.argsort(replica_indexes)]
    
    st1=np.zeros((n_snapshots,box_size**2,n_replicas),dtype=np.int32)
    n_replicas=0
    replica_indexes=np.array([],dtype=np.int32)
    for subdir, dirs, files in os.walk(path_to_array_replica_folder):
        for file in files:
            filepath = subdir + os.sep + file
            if "status1" in filepath:
                n_replicas += 1
                now_replica_index=int(filepath.split('/')[-2])
                replica_indexes= np.append(replica_indexes,now_replica_index)
                st1[:,:,n_replicas-1]=np.loadtxt(filepath)

    st1=st1[:,:,np.argsort(replica_indexes)]
    
    st2=np.zeros((n_snapshots,box_size**2,n_replicas),dtype=np.int32)
    n_replicas=0
    replica_indexes=np.array([],dtype=np.int32)
    for subdir, dirs, files in os.walk(path_to_array_replica_folder):
        for file in files:
            filepath = subdir + os.sep + file
            if "status2" in filepath:
                n_replicas += 1
                now_replica_index=int(filepath.split('/')[-2])
                replica_indexes= np.append(replica_indexes,now_replica_index)
                st2[:,:,n_replicas-1]=np.loadtxt(filepath)

    st2=st2[:,:,np.argsort(replica_indexes)]
    
    st3=np.zeros((n_snapshots,box_size**2,n_replicas),dtype=np.int32)
    n_replicas=0
    replica_indexes=np.array([],dtype=np.int32)
    for subdir, dirs, files in os.walk(path_to_array_replica_folder):
        for file in files:
            filepath = subdir + os.sep + file
            if "status3" in filepath:
                n_replicas += 1
                now_replica_index=int(filepath.split('/')[-2])
                replica_indexes= np.append(replica_indexes,now_replica_index)
                st3[:,:,n_replicas-1]=np.loadtxt(filepath)

    st3=st3[:,:,np.argsort(replica_indexes)]
    
    
    
    
    np.savez_compressed(path_to_array_replica_folder + 'analyzed_traj_loc',
        box_size=box_size,
        nativ_neigh_list=nativ_neigh_list,
        path_to_array_replica_folder=path_to_array_replica_folder,
        all_x_pos=all_x_pos,all_y_pos=all_y_pos,
        all_Nc=all_Nc,all_Ic=all_Ic,
        list_of_times=list_of_times,
        all_Nc_per_prot=all_Nc_per_prot,
        aa_present=aa_present,Nc_present=Nc_present,Ic_present=Ic_present,
        hw_present=hw_present,bw_present=bw_present,
        bulk_HB_present=bulk_HB_present,pho_HB_present=pho_HB_present,
        phi_HB_present=phi_HB_present,mix_HB_present=mix_HB_present,
        coop_bonds_present=coop_bonds_present,coop_bonds_PHO_present=coop_bonds_PHO_present,
        aa_counts=aa_counts,Nc_counts=Nc_counts,Ic_counts=Ic_counts,
        hw_counts=hw_counts,bw_counts=bw_counts,
        hb_b_counts=hb_b_counts,hb_pho_counts=hb_pho_counts,hb_phi_counts=hb_phi_counts,
        hb_mix_counts=hb_mix_counts,coop_b_counts=coop_b_counts,coop_pho_counts=coop_pho_counts,
        prot_counts=prot_counts,
        n_cells_tested=n_cells_tested,
        raw_data_folding=raw_data_folding,
        st0=st0,st1=st1,st2=st2,st3=st3)



def catenate_files_cluster_simple(npz_loaded_file_list):
    if type(npz_loaded_file_list) == list:
        # liq_cluster_label=npz_loaded_file_list[0]['liq_cluster_label']
        liq_cluster_size=npz_loaded_file_list[0]['liq_cluster_size']
        prot_belongs_to_liq_cluster=npz_loaded_file_list[0]['prot_belongs_to_liq_cluster']
        # sol_cluster_label=npz_loaded_file_list[0]['sol_cluster_label']
        sol_cluster_size=npz_loaded_file_list[0]['sol_cluster_size']
        prot_belongs_to_sol_cluster=npz_loaded_file_list[0]['prot_belongs_to_sol_cluster']
        # liq_cluster_label_nw=npz_loaded_file_list[0]['liq_cluster_label_nw']
        # liq_cluster_size_nw=npz_loaded_file_list[0]['liq_cluster_size_nw']
        # prot_belongs_to_liq_cluster_nw=npz_loaded_file_list[0]['prot_belongs_to_liq_cluster_nw']
        # sol_cluster_label_nw=npz_loaded_file_list[0]['sol_cluster_label_nw']
        # sol_cluster_size_nw=npz_loaded_file_list[0]['sol_cluster_size_nw']
        # prot_belongs_to_sol_cluster_nw=npz_loaded_file_list[0]['prot_belongs_to_sol_cluster_nw']
        tot_liq_cluster_sizes=npz_loaded_file_list[0]['tot_liq_cluster_sizes']
        tot_liq_cluster_volumes=npz_loaded_file_list[0]['tot_liq_cluster_volumes']
        tot_liq_cluster_volume_fractions=npz_loaded_file_list[0]['tot_liq_cluster_volume_fractions']
        tot_sol_cluster_sizes=npz_loaded_file_list[0]['tot_sol_cluster_sizes']
        tot_sol_cluster_volumes=npz_loaded_file_list[0]['tot_sol_cluster_volumes']
        tot_sol_cluster_volume_fractions=npz_loaded_file_list[0]['tot_sol_cluster_volume_fractions']
        
        cluster_hor_max_size=npz_loaded_file_list[0]['cluster_hor_max_size']
        cluster_ver_max_size=npz_loaded_file_list[0]['cluster_ver_max_size']
        sol_cluster_hor_max_size=npz_loaded_file_list[0]['sol_cluster_hor_max_size']
        sol_cluster_ver_max_size=npz_loaded_file_list[0]['sol_cluster_ver_max_size']
        liq_cluster_gyr_radius=npz_loaded_file_list[0]['liq_cluster_gyr_radius']
        sol_cluster_gyr_radius=npz_loaded_file_list[0]['sol_cluster_gyr_radius']
        for i_file in range(1,len(npz_loaded_file_list)):
            # liq_cluster_label=np.append(liq_cluster_label,npz_loaded_file_list[i_file]['liq_cluster_label'],axis=2)
            liq_cluster_size=np.append(liq_cluster_size,npz_loaded_file_list[i_file]['liq_cluster_size'],axis=1)
            prot_belongs_to_liq_cluster=np.append(prot_belongs_to_liq_cluster,npz_loaded_file_list[i_file]['prot_belongs_to_liq_cluster'],axis=1)
            # sol_cluster_label=np.append(sol_cluster_label,npz_loaded_file_list[i_file]['sol_cluster_label'],axis=2)
            sol_cluster_size=np.append(sol_cluster_size,npz_loaded_file_list[i_file]['sol_cluster_size'],axis=1)
            prot_belongs_to_sol_cluster=np.append(prot_belongs_to_sol_cluster,npz_loaded_file_list[i_file]['prot_belongs_to_sol_cluster'],axis=1)
            
            # liq_cluster_label_nw=np.append(liq_cluster_label_nw,npz_loaded_file_list[i_file]['liq_cluster_label_nw'],axis=2)
            # liq_cluster_size_nw=np.append(liq_cluster_size_nw,npz_loaded_file_list[i_file]['liq_cluster_size_nw'],axis=1)
            # prot_belongs_to_liq_cluster_nw=np.append(liq_cluster_size_nw,npz_loaded_file_list[i_file]['liq_cluster_size_nw'],axis=1)
            # sol_cluster_label_nw=np.append(liq_cluster_label_nw,npz_loaded_file_list[i_file]['liq_cluster_label_nw'],axis=2)
            # sol_cluster_size_nw=np.append(liq_cluster_size_nw,npz_loaded_file_list[i_file]['liq_cluster_size_nw'],axis=1)
            # prot_belongs_to_sol_cluster_nw=np.append(liq_cluster_size_nw,npz_loaded_file_list[i_fiaa_prevolumes'],axis=1)
            tot_liq_cluster_sizes=np.append(tot_liq_cluster_sizes,npz_loaded_file_list[i_file]['tot_liq_cluster_sizes'],axis=1)
            tot_liq_cluster_volumes=np.append(tot_liq_cluster_volumes,npz_loaded_file_list[i_file]['tot_liq_cluster_volumes'],axis=1)
            tot_liq_cluster_volume_fractions=np.append(tot_liq_cluster_volume_fractions,npz_loaded_file_list[i_file]['tot_liq_cluster_volume_fractions'],axis=1)
            tot_sol_cluster_sizes=np.append(tot_sol_cluster_sizes,npz_loaded_file_list[i_file]['tot_sol_cluster_sizes'],axis=1)
            tot_sol_cluster_volumes=np.append(tot_sol_cluster_volumes,npz_loaded_file_list[i_file]['tot_sol_cluster_volumes'],axis=1)
            tot_sol_cluster_volume_fractions=np.append(tot_sol_cluster_volume_fractions,npz_loaded_file_list[i_file]['tot_sol_cluster_volume_fractions'],axis=1)
            
            cluster_hor_max_size=np.append(cluster_hor_max_size,npz_loaded_file_list[i_file]['cluster_hor_max_size'],axis=0)
            cluster_ver_max_size=np.append(cluster_ver_max_size,npz_loaded_file_list[i_file]['cluster_ver_max_size'],axis=0)
            sol_cluster_hor_max_size=np.append(sol_cluster_hor_max_size,npz_loaded_file_list[i_file]['sol_cluster_hor_max_size'],axis=0)
            sol_cluster_ver_max_size=np.append(sol_cluster_ver_max_size,npz_loaded_file_list[i_file]['sol_cluster_ver_max_size'],axis=0)
            liq_cluster_gyr_radius=np.append(liq_cluster_gyr_radius,npz_loaded_file_list[i_file]['liq_cluster_gyr_radius'],axis=1)
            sol_cluster_gyr_radius=np.append(sol_cluster_gyr_radius,npz_loaded_file_list[i_file]['sol_cluster_gyr_radius'],axis=1)
    else:
        # liq_cluster_label=npz_loaded_file_list['liq_cluster_label']
        liq_cluster_size=npz_loaded_file_list['liq_cluster_size']
        prot_belongs_to_liq_cluster=npz_loaded_file_list['prot_belongs_to_liq_cluster']
        # sol_cluster_label=npz_loaded_file_list['sol_cluster_label']
        sol_cluster_size=npz_loaded_file_list['sol_cluster_size']
        prot_belongs_to_sol_cluster=npz_loaded_file_list['prot_belongs_to_sol_cluster']
        # liq_cluster_label_nw=npz_loaded_file_list['liq_cluster_label_nw']
        # liq_cluster_size_nw=npz_loaded_file_list['liq_cluster_size_nw']
        # prot_belongs_to_liq_cluster_nw=npz_loaded_file_list['prot_belongs_to_liq_cluster_nw']
        # sol_cluster_label_nw=npz_loaded_file_list['sol_cluster_label_nw']
        # sol_cluster_size_nw=npz_loaded_file_list['sol_cluster_size_nw']
        # prot_belongs_to_sol_cluster_nw=npz_loaded_file_list['prot_belongs_to_sol_cluster_nw']
        tot_liq_cluster_sizes=npz_loaded_file_list['tot_liq_cluster_sizes']
        tot_liq_cluster_volumes=npz_loaded_file_list['tot_liq_cluster_volumes']
        tot_liq_cluster_volume_fractions=npz_loaded_file_list['tot_liq_cluster_volume_fractions']
        tot_sol_cluster_sizes=npz_loaded_file_list['tot_sol_cluster_sizes']
        tot_sol_cluster_volumes=npz_loaded_file_list['tot_sol_cluster_volumes']
        tot_sol_cluster_volume_fractions=npz_loaded_file_list['tot_sol_cluster_volume_fractions']
        
        cluster_hor_max_size=npz_loaded_file_list['cluster_hor_max_size']
        cluster_ver_max_size=npz_loaded_file_list['cluster_ver_max_size']
        liq_cluster_gyr_radius=npz_loaded_file_list['liq_cluster_gyr_radius']
        sol_cluster_gyr_radius=npz_loaded_file_list['sol_cluster_gyr_radius']
        sol_cluster_hor_max_size=npz_loaded_file_list['sol_cluster_hor_max_size']
        sol_cluster_ver_max_size=npz_loaded_file_list['sol_cluster_ver_max_size']
        
    return liq_cluster_size,prot_belongs_to_liq_cluster,\
        sol_cluster_size,prot_belongs_to_sol_cluster,\
        tot_liq_cluster_sizes,tot_liq_cluster_volumes,tot_liq_cluster_volume_fractions,\
        tot_sol_cluster_sizes,tot_sol_cluster_volumes,tot_sol_cluster_volume_fractions,\
        cluster_hor_max_size,cluster_ver_max_size,\
        sol_cluster_hor_max_size,sol_cluster_ver_max_size,\
        liq_cluster_gyr_radius,sol_cluster_gyr_radius
        


def catenate_files_cluster(npz_loaded_file_list):
    if type(npz_loaded_file_list) == list:
        liq_cluster_label=npz_loaded_file_list[0]['liq_cluster_label']
        liq_cluster_size=npz_loaded_file_list[0]['liq_cluster_size']
        prot_belongs_to_liq_cluster=npz_loaded_file_list[0]['prot_belongs_to_liq_cluster']
        sol_cluster_label=npz_loaded_file_list[0]['sol_cluster_label']
        sol_cluster_size=npz_loaded_file_list[0]['sol_cluster_size']
        prot_belongs_to_sol_cluster=npz_loaded_file_list[0]['prot_belongs_to_sol_cluster']
        liq_cluster_label_nw=npz_loaded_file_list[0]['liq_cluster_label_nw']
        liq_cluster_size_nw=npz_loaded_file_list[0]['liq_cluster_size_nw']
        prot_belongs_to_liq_cluster_nw=npz_loaded_file_list[0]['prot_belongs_to_liq_cluster_nw']
        sol_cluster_label_nw=npz_loaded_file_list[0]['sol_cluster_label_nw']
        sol_cluster_size_nw=npz_loaded_file_list[0]['sol_cluster_size_nw']
        prot_belongs_to_sol_cluster_nw=npz_loaded_file_list[0]['prot_belongs_to_sol_cluster_nw']
        tot_liq_cluster_sizes=npz_loaded_file_list[0]['tot_liq_cluster_sizes']
        tot_liq_cluster_volumes=npz_loaded_file_list[0]['tot_liq_cluster_volumes']
        tot_liq_cluster_volume_fractions=npz_loaded_file_list[0]['tot_liq_cluster_volume_fractions']
        tot_sol_cluster_sizes=npz_loaded_file_list[0]['tot_sol_cluster_sizes']
        tot_sol_cluster_volumes=npz_loaded_file_list[0]['tot_sol_cluster_volumes']
        tot_sol_cluster_volume_fractions=npz_loaded_file_list[0]['tot_sol_cluster_volume_fractions']
        cluster_hor_max_size=npz_loaded_file_list[0]['cluster_hor_max_size']
        cluster_ver_max_size=npz_loaded_file_list[0]['cluster_ver_max_size']
        sol_cluster_hor_max_size=npz_loaded_file_list[0]['sol_cluster_hor_max_size']
        sol_cluster_ver_max_size=npz_loaded_file_list[0]['sol_cluster_ver_max_size']
        liq_cluster_gyr_radius=npz_loaded_file_list[0]['liq_cluster_gyr_radius']
        sol_cluster_gyr_radius=npz_loaded_file_list[0]['sol_cluster_gyr_radius']
        for i_file in range(1,len(npz_loaded_file_list)):
            liq_cluster_label=np.append(liq_cluster_label,npz_loaded_file_list[i_file]['liq_cluster_label'],axis=2)
            liq_cluster_size=np.append(liq_cluster_size,npz_loaded_file_list[i_file]['liq_cluster_size'],axis=1)
            prot_belongs_to_liq_cluster=np.append(prot_belongs_to_liq_cluster,npz_loaded_file_list[i_file]['prot_belongs_to_liq_cluster'],axis=1)
            sol_cluster_label=np.append(sol_cluster_label,npz_loaded_file_list[i_file]['sol_cluster_label'],axis=2)
            sol_cluster_size=np.append(sol_cluster_size,npz_loaded_file_list[i_file]['sol_cluster_size'],axis=1)
            prot_belongs_to_sol_cluster=np.append(prot_belongs_to_sol_cluster,npz_loaded_file_list[i_file]['prot_belongs_to_sol_cluster'],axis=1)
            
            liq_cluster_label_nw=np.append(liq_cluster_label_nw,npz_loaded_file_list[i_file]['liq_cluster_label_nw'],axis=2)
            liq_cluster_size_nw=np.append(liq_cluster_size_nw,npz_loaded_file_list[i_file]['liq_cluster_size_nw'],axis=1)
            prot_belongs_to_liq_cluster_nw=np.append(liq_cluster_size_nw,npz_loaded_file_list[i_file]['liq_cluster_size_nw'],axis=1)
            sol_cluster_label_nw=np.append(liq_cluster_label_nw,npz_loaded_file_list[i_file]['liq_cluster_label_nw'],axis=2)
            sol_cluster_size_nw=np.append(liq_cluster_size_nw,npz_loaded_file_list[i_file]['liq_cluster_size_nw'],axis=1)
            prot_belongs_to_sol_cluster_nw=np.append(liq_cluster_size_nw,npz_loaded_file_list[i_file]['liq_cluster_size_nw'],axis=1)
            
            tot_liq_cluster_sizes=np.append(tot_liq_cluster_sizes,npz_loaded_file_list[i_file]['tot_liq_cluster_sizes'],axis=1)
            tot_liq_cluster_volumes=np.append(tot_liq_cluster_volumes,npz_loaded_file_list[i_file]['tot_liq_cluster_volumes'],axis=1)
            tot_liq_cluster_volume_fractions=np.append(tot_liq_cluster_volume_fractions,npz_loaded_file_list[i_file]['tot_liq_cluster_volume_fractions'],axis=1)
            tot_sol_cluster_sizes=np.append(tot_sol_cluster_sizes,npz_loaded_file_list[i_file]['tot_sol_cluster_sizes'],axis=1)
            tot_sol_cluster_volumes=np.append(tot_sol_cluster_volumes,npz_loaded_file_list[i_file]['tot_sol_cluster_volumes'],axis=1)
            tot_sol_cluster_volume_fractions=np.append(tot_sol_cluster_volume_fractions,npz_loaded_file_list[i_file]['tot_sol_cluster_volume_fractions'],axis=1)
            
            cluster_hor_max_size=np.append(cluster_hor_max_size,npz_loaded_file_list[i_file]['cluster_hor_max_size'],axis=0)
            cluster_ver_max_size=np.append(cluster_ver_max_size,npz_loaded_file_list[i_file]['cluster_ver_max_size'],axis=0)
            sol_cluster_hor_max_size=np.append(sol_cluster_hor_max_size,npz_loaded_file_list[i_file]['sol_cluster_hor_max_size'],axis=0)
            sol_cluster_ver_max_size=np.append(sol_cluster_ver_max_size,npz_loaded_file_list[i_file]['sol_cluster_ver_max_size'],axis=0)
            liq_cluster_gyr_radius=np.append(liq_cluster_gyr_radius,npz_loaded_file_list[i_file]['liq_cluster_gyr_radius'],axis=1)
            sol_cluster_gyr_radius=np.append(sol_cluster_gyr_radius,npz_loaded_file_list[i_file]['sol_cluster_gyr_radius'],axis=1)

    else:
        liq_cluster_label=npz_loaded_file_list['liq_cluster_label']
        liq_cluster_size=npz_loaded_file_list['liq_cluster_size']
        prot_belongs_to_liq_cluster=npz_loaded_file_list['prot_belongs_to_liq_cluster']
        sol_cluster_label=npz_loaded_file_list['sol_cluster_label']
        sol_cluster_size=npz_loaded_file_list['sol_cluster_size']
        prot_belongs_to_sol_cluster=npz_loaded_file_list['prot_belongs_to_sol_cluster']
        liq_cluster_label_nw=npz_loaded_file_list['liq_cluster_label_nw']
        liq_cluster_size_nw=npz_loaded_file_list['liq_cluster_size_nw']
        prot_belongs_to_liq_cluster_nw=npz_loaded_file_list['prot_belongs_to_liq_cluster_nw']
        sol_cluster_label_nw=npz_loaded_file_list['sol_cluster_label_nw']
        sol_cluster_size_nw=npz_loaded_file_list['sol_cluster_size_nw']
        prot_belongs_to_sol_cluster_nw=npz_loaded_file_list['prot_belongs_to_sol_cluster_nw']
        tot_liq_cluster_sizes=npz_loaded_file_list['tot_liq_cluster_sizes']
        tot_liq_cluster_volumes=npz_loaded_file_list['tot_liq_cluster_volumes']
        tot_liq_cluster_volume_fractions=npz_loaded_file_list['tot_liq_cluster_volume_fractions']
        tot_sol_cluster_sizes=npz_loaded_file_list['tot_sol_cluster_sizes']
        tot_sol_cluster_volumes=npz_loaded_file_list['tot_sol_cluster_volumes']
        tot_sol_cluster_volume_fractions=npz_loaded_file_list['tot_sol_cluster_volume_fractions']
        cluster_hor_max_size=npz_loaded_file_list['cluster_hor_max_size']
        cluster_ver_max_size=npz_loaded_file_list['cluster_ver_max_size']
        liq_cluster_gyr_radius=npz_loaded_file_list['liq_cluster_gyr_radius']
        sol_cluster_gyr_radius=npz_loaded_file_list['sol_cluster_gyr_radius']
        sol_cluster_hor_max_size=npz_loaded_file_list['sol_cluster_hor_max_size']
        sol_cluster_ver_max_size=npz_loaded_file_list['sol_cluster_ver_max_size']
        
        
        
    return liq_cluster_label,liq_cluster_size,prot_belongs_to_liq_cluster,\
        sol_cluster_label,sol_cluster_size,prot_belongs_to_sol_cluster,\
        liq_cluster_label_nw,liq_cluster_size_nw,prot_belongs_to_liq_cluster_nw,\
        sol_cluster_label_nw,sol_cluster_size_nw,prot_belongs_to_sol_cluster_nw,\
        tot_liq_cluster_sizes,tot_liq_cluster_volumes,tot_liq_cluster_volume_fractions,\
        tot_sol_cluster_sizes,tot_sol_cluster_volumes,tot_sol_cluster_volume_fractions,\
        cluster_hor_max_size,cluster_ver_max_size,\
        sol_cluster_hor_max_size,sol_cluster_ver_max_size,\
        liq_cluster_gyr_radius,sol_cluster_gyr_radius


def catenate_files_local_all_savings(npz_loaded_file_list):
    if type(npz_loaded_file_list) == list:
        box_size=int(npz_loaded_file_list[0]['box_size'])
        nativ_neigh_list=npz_loaded_file_list[0]['nativ_neigh_list']
        path_to_array_replica_folder=npz_loaded_file_list[0]['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list[0]['all_x_pos']
        all_y_pos=npz_loaded_file_list[0]['all_y_pos']
        all_Nc=npz_loaded_file_list[0]['all_Nc']
        all_Ic=npz_loaded_file_list[0]['all_Ic']
        list_of_times=npz_loaded_file_list[0]['list_of_times']
        all_Nc_per_prot=npz_loaded_file_list[0]['all_Nc_per_prot']
        
        aa_present=npz_loaded_file_list[0]['aa_present']
        Nc_present=npz_loaded_file_list[0]['Nc_present']
        Ic_present=npz_loaded_file_list[0]['Ic_present']
        hw_present=npz_loaded_file_list[0]['hw_present']
        bw_present=npz_loaded_file_list[0]['bw_present']
        bulk_HB_present=npz_loaded_file_list[0]['bulk_HB_present']
        pho_HB_present=npz_loaded_file_list[0]['pho_HB_present']
        phi_HB_present=npz_loaded_file_list[0]['phi_HB_present']
        mix_HB_present=npz_loaded_file_list[0]['mix_HB_present']
        coop_bonds_present=npz_loaded_file_list[0]['coop_bonds_present']
        coop_bonds_PHO_present=npz_loaded_file_list[0]['coop_bonds_PHO_present']
        
        aa_counts=npz_loaded_file_list[0]['aa_counts']
        Nc_counts=npz_loaded_file_list[0]['Nc_counts']
        Ic_counts=npz_loaded_file_list[0]['Ic_counts']
        hw_counts=npz_loaded_file_list[0]['hw_counts']
        bw_counts=npz_loaded_file_list[0]['bw_counts']
        hb_b_counts=npz_loaded_file_list[0]['hb_b_counts']
        hb_pho_counts=npz_loaded_file_list[0]['hb_pho_counts']
        hb_phi_counts=npz_loaded_file_list[0]['hb_phi_counts']
        hb_mix_counts=npz_loaded_file_list[0]['hb_mix_counts']
        coop_b_counts=npz_loaded_file_list[0]['coop_b_counts']
        coop_pho_counts=npz_loaded_file_list[0]['coop_pho_counts']
        
        prot_counts=npz_loaded_file_list[0]['prot_counts']
        
        n_cells_tested=npz_loaded_file_list[0]['n_cells_tested']
        raw_data_folding=npz_loaded_file_list[0]['raw_data_folding']
        
        st0=npz_loaded_file_list[0]['st0']
        st1=npz_loaded_file_list[0]['st1']
        st2=npz_loaded_file_list[0]['st2']
        st3=npz_loaded_file_list[0]['st3']
        
        for i_file in range(1,len(npz_loaded_file_list)):
            all_x_pos=np.append(all_x_pos,npz_loaded_file_list[i_file]['all_x_pos'],axis=2)
            all_y_pos=np.append(all_y_pos,npz_loaded_file_list[i_file]['all_y_pos'],axis=2)
            all_Nc=np.append(all_Nc,npz_loaded_file_list[i_file]['all_Nc'],axis=1)
            all_Ic=np.append(all_Ic,npz_loaded_file_list[i_file]['all_Ic'],axis=1)
            list_of_times=np.append(list_of_times,npz_loaded_file_list[i_file]['list_of_times']+list_of_times[-1])
            all_Nc_per_prot=np.append(all_Nc_per_prot,npz_loaded_file_list[i_file]['all_Nc_per_prot'],axis=2)

            aa_present=np.append(aa_present,npz_loaded_file_list[i_file]['aa_present'],axis=2)
            Nc_present=np.append(Nc_present,npz_loaded_file_list[i_file]['Nc_present'],axis=2)
            Ic_present=np.append(Ic_present,npz_loaded_file_list[i_file]['Ic_present'],axis=2)
            hw_present=np.append(hw_present,npz_loaded_file_list[i_file]['hw_present'],axis=2)
            bw_present=np.append(bw_present,npz_loaded_file_list[i_file]['bw_present'],axis=2)
            bulk_HB_present=np.append(bulk_HB_present,npz_loaded_file_list[i_file]['bulk_HB_present'],axis=2)
            pho_HB_present=np.append(pho_HB_present,npz_loaded_file_list[i_file]['pho_HB_present'],axis=2)
            phi_HB_present=np.append(phi_HB_present,npz_loaded_file_list[i_file]['phi_HB_present'],axis=2)
            mix_HB_present=np.append(mix_HB_present,npz_loaded_file_list[i_file]['mix_HB_present'],axis=2)
            coop_bonds_present=np.append(coop_bonds_present,npz_loaded_file_list[i_file]['coop_bonds_present'],axis=2)
            coop_bonds_PHO_present=np.append(coop_bonds_PHO_present,npz_loaded_file_list[i_file]['coop_bonds_PHO_present'],axis=2)
            
            aa_counts=np.append(aa_counts,npz_loaded_file_list[i_file]['aa_counts'],axis=2)
            Nc_counts=np.append(Nc_counts,npz_loaded_file_list[i_file]['Nc_counts'],axis=2)
            Ic_counts=np.append(Ic_counts,npz_loaded_file_list[i_file]['Ic_counts'],axis=2)
            hw_counts=np.append(hw_counts,npz_loaded_file_list[i_file]['hw_counts'],axis=2)
            bw_counts=np.append(bw_counts,npz_loaded_file_list[i_file]['bw_counts'],axis=2)
            hb_b_counts=np.append(hb_b_counts,npz_loaded_file_list[i_file]['hb_b_counts'],axis=2)
            hb_pho_counts=np.append(hb_pho_counts,npz_loaded_file_list[i_file]['hb_pho_counts'],axis=2)
            hb_phi_counts=np.append(hb_phi_counts,npz_loaded_file_list[i_file]['hb_phi_counts'],axis=2)
            hb_mix_counts=np.append(hb_mix_counts,npz_loaded_file_list[i_file]['hb_mix_counts'],axis=2)
            coop_b_counts=np.append(coop_b_counts,npz_loaded_file_list[i_file]['coop_b_counts'],axis=2)
            coop_pho_counts=np.append(coop_pho_counts,npz_loaded_file_list[i_file]['coop_pho_counts'],axis=2)
            
            prot_counts=np.append(prot_counts,npz_loaded_file_list[i_file]['prot_counts'],axis=2)
            
            raw_data_folding=np.append(raw_data_folding,npz_loaded_file_list[i_file]['raw_data_folding'],axis=0)
            
            st0=np.append(st0,npz_loaded_file_list[i_file]['st0'],axis=0)
            st1=np.append(st1,npz_loaded_file_list[i_file]['st1'],axis=0)
            st2=np.append(st2,npz_loaded_file_list[i_file]['st2'],axis=0)
            st3=np.append(st3,npz_loaded_file_list[i_file]['st3'],axis=0)
    else:
        box_size=int(npz_loaded_file_list['box_size'])
        nativ_neigh_list=npz_loaded_file_list['nativ_neigh_list']
        path_to_array_replica_folder=npz_loaded_file_list['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list['all_x_pos']
        all_y_pos=npz_loaded_file_list['all_y_pos']
        all_Nc=npz_loaded_file_list['all_Nc']
        all_Ic=npz_loaded_file_list['all_Ic']
        list_of_times=npz_loaded_file_list['list_of_times']
        all_Nc_per_prot=npz_loaded_file_list['all_Nc_per_prot']
        aa_present=npz_loaded_file_list['aa_present']
        Nc_present=npz_loaded_file_list['Nc_present']
        Ic_present=npz_loaded_file_list['Ic_present']
        hw_present=npz_loaded_file_list['hw_present']
        bw_present=npz_loaded_file_list['bw_present']
        bulk_HB_present=npz_loaded_file_list['bulk_HB_present']
        pho_HB_present=npz_loaded_file_list['pho_HB_present']
        phi_HB_present=npz_loaded_file_list['phi_HB_present']
        mix_HB_present=npz_loaded_file_list['mix_HB_present']
        coop_bonds_present=npz_loaded_file_list['coop_bonds_present']
        coop_bonds_PHO_present=npz_loaded_file_list['coop_bonds_PHO_present']
        
        aa_counts=npz_loaded_file_list['aa_counts']
        Nc_counts=npz_loaded_file_list['Nc_counts']
        Ic_counts=npz_loaded_file_list['Ic_counts']
        hw_counts=npz_loaded_file_list['hw_counts']
        bw_counts=npz_loaded_file_list['bw_counts']
        hb_b_counts=npz_loaded_file_list['hb_b_counts']
        hb_pho_counts=npz_loaded_file_list['hb_pho_counts']
        hb_phi_counts=npz_loaded_file_list['hb_phi_counts']
        hb_mix_counts=npz_loaded_file_list['hb_mix_counts']
        coop_b_counts=npz_loaded_file_list['coop_b_counts']
        coop_pho_counts=npz_loaded_file_list['coop_pho_counts']
        
        prot_counts=npz_loaded_file_list['prot_counts']
        
        n_cells_tested=npz_loaded_file_list['n_cells_tested']
        raw_data_folding=npz_loaded_file_list['raw_data_folding']
        
        st0=npz_loaded_file_list['st0']
        st1=npz_loaded_file_list['st1']
        st2=npz_loaded_file_list['st2']
        st3=npz_loaded_file_list['st3']
        
    return box_size,nativ_neigh_list,path_to_array_replica_folder,\
    all_x_pos,all_y_pos,all_Nc,all_Ic,list_of_times,all_Nc_per_prot,\
    aa_present,Nc_present,Ic_present,hw_present,bw_present,\
    bulk_HB_present,pho_HB_present,phi_HB_present,mix_HB_present,coop_bonds_present,coop_bonds_PHO_present,\
    aa_counts,Nc_counts,Ic_counts,hw_counts,bw_counts,\
    hb_b_counts,hb_pho_counts,hb_phi_counts,hb_mix_counts,coop_b_counts,coop_pho_counts,\
    prot_counts,\
    n_cells_tested,\
    raw_data_folding,\
    st0,st1,st2,st3



def catenate_files_local_simple(npz_loaded_file_list):
    if type(npz_loaded_file_list) == list:
        box_size=int(npz_loaded_file_list[0]['box_size'])
        nativ_neigh_list=npz_loaded_file_list[0]['nativ_neigh_list']
        path_to_array_replica_folder=npz_loaded_file_list[0]['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list[0]['all_x_pos']
        all_y_pos=npz_loaded_file_list[0]['all_y_pos']
        all_Nc=npz_loaded_file_list[0]['all_Nc']
        all_Ic=npz_loaded_file_list[0]['all_Ic']
        list_of_times=npz_loaded_file_list[0]['list_of_times']
                  
        raw_data_folding=npz_loaded_file_list[0]['raw_data_folding']
        
        # st0=npz_loaded_file_list[0]['st0']
        # st1=npz_loaded_file_list[0]['st1']
        # st2=npz_loaded_file_list[0]['st2']
        # st3=npz_loaded_file_list[0]['st3']
        
        for i_file in range(1,len(npz_loaded_file_list)):
            all_x_pos=np.append(all_x_pos,npz_loaded_file_list[i_file]['all_x_pos'],axis=2)
            all_y_pos=np.append(all_y_pos,npz_loaded_file_list[i_file]['all_y_pos'],axis=2)
            all_Nc=np.append(all_Nc,npz_loaded_file_list[i_file]['all_Nc'],axis=1)
            all_Ic=np.append(all_Ic,npz_loaded_file_list[i_file]['all_Ic'],axis=1)
            list_of_times=np.append(list_of_times,npz_loaded_file_list[i_file]['list_of_times']+list_of_times[-1])
                               
            raw_data_folding=np.append(raw_data_folding,npz_loaded_file_list[i_file]['raw_data_folding'],axis=0)
            
            # st0=np.append(st0,npz_loaded_file_list[i_file]['st0'],axis=0)
            # st1=np.append(st1,npz_loaded_file_list[i_file]['st1'],axis=0)
            # st2=np.append(st2,npz_loaded_file_list[i_file]['st2'],axis=0)
            # st3=np.append(st3,npz_loaded_file_list[i_file]['st3'],axis=0)
    else:
        box_size=int(npz_loaded_file_list['box_size'])
        nativ_neigh_list=npz_loaded_file_list['nativ_neigh_list']
        path_to_array_replica_folder=npz_loaded_file_list['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list['all_x_pos']
        all_y_pos=npz_loaded_file_list['all_y_pos']
        all_Nc=npz_loaded_file_list['all_Nc']
        all_Ic=npz_loaded_file_list['all_Ic']
        list_of_times=npz_loaded_file_list['list_of_times']

                
        raw_data_folding=npz_loaded_file_list['raw_data_folding']
        
        # st0=npz_loaded_file_list['st0']
        # st1=npz_loaded_file_list['st1']
        # st2=npz_loaded_file_list['st2']
        # st3=npz_loaded_file_list['st3']
        
    return box_size,nativ_neigh_list,path_to_array_replica_folder,\
    all_x_pos,all_y_pos,all_Nc,all_Ic,list_of_times,\
    raw_data_folding#,\
    # st0,st1,st2,st3
    
    




    


def get_average_in_t_window(disc_phi,n_cells_tested,\
    aa_counts,Nc_counts,Ic_counts,hw_counts,bw_counts,
    in_t,fi_t,n_native_bonds_per_protein,n_monomers):
    disc_aa_counts=int(n_cells_tested*disc_phi)
    n_native_bonds_per_monomer=n_native_bonds_per_protein/n_monomers
    
    now_mask_aa_counts_g_disc=aa_counts>disc_aa_counts
    now_mask_aa_counts_le_disc=aa_counts<=disc_aa_counts
    now_mask_aa_counts_all=np.ones(now_mask_aa_counts_g_disc.shape,dtype=np.bool_)
    now_mask_aa_counts_g_disc[:,:,:in_t,:]=False
    now_mask_aa_counts_g_disc[:,:,fi_t:,:]=False
    now_mask_aa_counts_le_disc[:,:,:in_t,:]=False
    now_mask_aa_counts_le_disc[:,:,fi_t:,:]=False
    now_mask_aa_counts_all[:,:,:in_t,:]=False
    now_mask_aa_counts_all[:,:,fi_t:,:]=False
    
    system_fracs=[now_mask_aa_counts_g_disc[:,:,in_t:fi_t,:].mean(axis=(0,1)),#/(n_replicas*box_size**2),
                  now_mask_aa_counts_le_disc[:,:,in_t:fi_t,:].mean(axis=(0,1))]

    ave_system_fracs=[system_fracs[0].mean(),system_fracs[1].mean()]
    ave_system_fracs_std=[system_fracs[0].std(),system_fracs[1].std()]
    
    conc=[aa_counts[now_mask_aa_counts_all]/n_cells_tested,
        aa_counts[now_mask_aa_counts_g_disc]/n_cells_tested,
        aa_counts[now_mask_aa_counts_le_disc]/n_cells_tested]
    
    conc_hw=[hw_counts[now_mask_aa_counts_all]/n_cells_tested,
        hw_counts[now_mask_aa_counts_g_disc]/n_cells_tested,
        hw_counts[now_mask_aa_counts_le_disc]/n_cells_tested]
    
    conc_bw=[bw_counts[now_mask_aa_counts_all]/n_cells_tested,
        bw_counts[now_mask_aa_counts_g_disc]/n_cells_tested,
        bw_counts[now_mask_aa_counts_le_disc]/n_cells_tested]
    
    
    ave_conc=[
        aa_counts[now_mask_aa_counts_all].mean()/n_cells_tested,
        aa_counts[now_mask_aa_counts_g_disc].mean()/n_cells_tested,
        aa_counts[now_mask_aa_counts_le_disc].mean()/n_cells_tested]
    ave_conc_std=[
        aa_counts[now_mask_aa_counts_all].std()/n_cells_tested,
        aa_counts[now_mask_aa_counts_g_disc].std()/n_cells_tested,
        aa_counts[now_mask_aa_counts_le_disc].std()/n_cells_tested]
    
    ave_hw_conc=[
        hw_counts[now_mask_aa_counts_all].mean()/n_cells_tested,
        hw_counts[now_mask_aa_counts_g_disc].mean()/n_cells_tested,
        hw_counts[now_mask_aa_counts_le_disc].mean()/n_cells_tested]
    ave_hw_conc_std=[
        hw_counts[now_mask_aa_counts_all].std()/n_cells_tested,
        hw_counts[now_mask_aa_counts_g_disc].std()/n_cells_tested,
        hw_counts[now_mask_aa_counts_le_disc].std()/n_cells_tested]
    
    ave_bw_conc=[
        bw_counts[now_mask_aa_counts_all].mean()/n_cells_tested,
        bw_counts[now_mask_aa_counts_g_disc].mean()/n_cells_tested,
        bw_counts[now_mask_aa_counts_le_disc].mean()/n_cells_tested]
    ave_bw_conc_std=[
        bw_counts[now_mask_aa_counts_all].std()/n_cells_tested,
        bw_counts[now_mask_aa_counts_g_disc].std()/n_cells_tested,
        bw_counts[now_mask_aa_counts_le_disc].std()/n_cells_tested]
    
   
    Nc_degree=[
        np.divide(Nc_counts[now_mask_aa_counts_all][aa_counts[now_mask_aa_counts_all]>0],
            aa_counts[now_mask_aa_counts_all][aa_counts[now_mask_aa_counts_all]>0])/(2*n_native_bonds_per_monomer),
        np.divide(Nc_counts[now_mask_aa_counts_g_disc][aa_counts[now_mask_aa_counts_g_disc]>0],
            aa_counts[now_mask_aa_counts_g_disc][aa_counts[now_mask_aa_counts_g_disc]>0])/(2*n_native_bonds_per_monomer),
        np.divide(Nc_counts[now_mask_aa_counts_le_disc][aa_counts[now_mask_aa_counts_le_disc]>0],
            aa_counts[now_mask_aa_counts_le_disc][aa_counts[now_mask_aa_counts_le_disc]>0])/(2*n_native_bonds_per_monomer)]
    
    ave_Nc_degree=[Nc_degree[0].mean(),Nc_degree[1].mean(),Nc_degree[2].mean()]
    ave_Nc_degree_std=[Nc_degree[0].std(),Nc_degree[1].std(),Nc_degree[2].std()]
    
    Ic_degree=[
        np.divide(Ic_counts[now_mask_aa_counts_all][aa_counts[now_mask_aa_counts_all]>0],
            aa_counts[now_mask_aa_counts_all][aa_counts[now_mask_aa_counts_all]>0])/(2*n_native_bonds_per_monomer),
        np.divide(Ic_counts[now_mask_aa_counts_g_disc][aa_counts[now_mask_aa_counts_g_disc]>0],
            aa_counts[now_mask_aa_counts_g_disc][aa_counts[now_mask_aa_counts_all]>0])/(2*n_native_bonds_per_monomer),
        np.divide(Ic_counts[now_mask_aa_counts_le_disc][aa_counts[now_mask_aa_counts_le_disc]>0],
            aa_counts[now_mask_aa_counts_le_disc][aa_counts[now_mask_aa_counts_le_disc]>0])/(2*n_native_bonds_per_monomer)]
    
    ave_Ic_degree=[Ic_degree[0].mean(),Ic_degree[1].mean(),Ic_degree[2].mean()]
    ave_Ic_degree_std=[Ic_degree[0].std(),Ic_degree[1].std(),Ic_degree[2].std()]
    
    # aux_Nc_counts=Nc_counts
    # aux_aa_counts=aa_counts
    # aux_Nc_counts[aa_counts==0]=0
    # aux_aa_counts[aa_counts==0]=1
    # Nc_degree=[
    #     np.divide(aux_Nc_counts[now_mask_aa_counts_all],aux_aa_counts[now_mask_aa_counts_all])/(2*n_native_bonds_per_monomer),
    #     np.divide(aux_Nc_counts[now_mask_aa_counts_g_disc],aux_aa_counts[now_mask_aa_counts_g_disc])/(2*n_native_bonds_per_monomer),
    #     np.divide(aux_Nc_counts[now_mask_aa_counts_le_disc],aux_aa_counts[now_mask_aa_counts_le_disc])/(2*n_native_bonds_per_monomer),]
    



    return ave_system_fracs,ave_system_fracs_std,ave_conc,ave_conc_std,\
        ave_Nc_degree,ave_Nc_degree_std,ave_Ic_degree,ave_Ic_degree_std,\
            ave_hw_conc,ave_hw_conc_std,ave_bw_conc,ave_bw_conc_std,\
            conc,conc_hw,conc_bw,Nc_degree,Ic_degree








@njit
def convert_HB_to_mol_states(all_hor_HB,all_ver_HB,all_coop_bonds,all_coop_bonds_PHO,box_size):
    n_snapshots=all_hor_HB.shape[0]
    n_replicas=all_hor_HB.shape[2]
    rx=np.zeros((box_size**2),dtype=np.int32)
    ry=np.zeros((box_size**2),dtype=np.int32)
    i_site=0
    for ix in range(box_size):
        for iy in range(box_size):
            rx[i_site]=ix+1
            ry[i_site]=iy+1
            i_site+=1
            
            
    bulk_HB_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas),dtype=np.int32)
    pho_HB_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas),dtype=np.int32)
    phi_HB_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas),dtype=np.int32)
    mix_HB_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas),dtype=np.int32)
    coop_bonds_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas),dtype=np.int32)
    coop_bonds_PHO_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas),dtype=np.int32)
    
    
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            hor_HB=all_hor_HB[i_t,:,i_replica]
            ver_HB=all_ver_HB[i_t,:,i_replica]
            coop_bonds=all_coop_bonds[i_t,:,i_replica]
            coop_bonds_PHO=all_coop_bonds_PHO[i_t,:,i_replica]
            for i_site in range(box_size**2):
                now_rx=rx[i_site]
                now_ry=ry[i_site]
                now_now_rxp1=((now_rx+1+box_size-1)%box_size)+1
                now_now_rym1=((now_ry-1+box_size-1)%box_size)+1
                # check for bulk
                if (hor_HB[i_site]==1):
                    bulk_HB_present[now_rx,now_ry,i_t,i_replica]+=1
                    bulk_HB_present[now_now_rxp1,now_ry,i_t,i_replica]+=1
                if (ver_HB[i_site]==1):
                    bulk_HB_present[now_rx,now_ry,i_t,i_replica]+=1
                    bulk_HB_present[now_rx,now_now_rym1,i_t,i_replica]+=1
                # check for PHO (2)
                if (hor_HB[i_site]==2):
                    pho_HB_present[now_rx,now_ry,i_t,i_replica]+=1
                    pho_HB_present[now_now_rxp1,now_ry,i_t,i_replica]+=1
                if (ver_HB[i_site]==2):
                    pho_HB_present[now_rx,now_ry,i_t,i_replica]+=1
                    pho_HB_present[now_rx,now_now_rym1,i_t,i_replica]+=1
                # check for PHI (3)
                if (hor_HB[i_site]==3):
                    phi_HB_present[now_rx,now_ry,i_t,i_replica]+=1
                    phi_HB_present[now_now_rxp1,now_ry,i_t,i_replica]+=1
                if (ver_HB[i_site]==3):
                    phi_HB_present[now_rx,now_ry,i_t,i_replica]+=1
                    phi_HB_present[now_rx,now_now_rym1,i_t,i_replica]+=1
                # check for MIX (3)
                if (hor_HB[i_site]==4):
                    mix_HB_present[now_rx,now_ry,i_t,i_replica]+=1
                    mix_HB_present[now_now_rxp1,now_ry,i_t,i_replica]+=1
                if (ver_HB[i_site]==4):
                    mix_HB_present[now_rx,now_ry,i_t,i_replica]+=1
                    mix_HB_present[now_rx,now_now_rym1,i_t,i_replica]+=1
                if (coop_bonds[i_site]>0):
                    coop_bonds_present[now_rx,now_ry,i_t,i_replica]+=coop_bonds[i_site]
                if (coop_bonds_PHO[i_site]>0):
                    coop_bonds_PHO_present[now_rx,now_ry,i_t,i_replica]+=coop_bonds_PHO[i_site]

    return bulk_HB_present,pho_HB_present,phi_HB_present,mix_HB_present,coop_bonds_present,coop_bonds_PHO_present





@njit
def calculate_all_local(all_x_pos,all_y_pos,box_size,test_radius,nativ_neigh_list, \
    bulk_HB_present,pho_HB_present,phi_HB_present,mix_HB_present,coop_bonds_present,coop_bonds_PHO_present):
    n_monomers=all_x_pos.shape[0]
    n_proteins=all_x_pos.shape[1]
    n_snapshots=all_x_pos.shape[2]
    n_replicas=all_x_pos.shape[3]
    
    # count number of cells involved in calcualtion for given test_radius
    n_cells_tested=0
    now_center_x=1
    now_center_y=1
    for i_x in range(-test_radius,test_radius+1):
        for i_y in range(-test_radius,test_radius+1):
            now_x_ind=((now_center_x+i_x+box_size-1)%box_size)+1
            now_y_ind=((now_center_y+i_y+box_size-1)%box_size)+1
            dx = box_size//2 - abs(abs(now_x_ind-now_center_x) - box_size//2)
            dy = box_size//2 - abs(abs(now_y_ind-now_center_y) - box_size//2)
            if (dx**2+dy**2 <= test_radius**2):
                n_cells_tested+=1
    # count how many cells inside radius are filled
    # repeat for every replica and every snapshot
    #
    aa_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    
    aa_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas), dtype=np.bool_)
    Nc_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas), dtype=np.int32)
    Ic_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas), dtype=np.int32)
    hw_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas), dtype=np.bool_)
    bw_present=np.zeros((box_size+1,box_size+1,n_snapshots,n_replicas), dtype=np.bool_)
    
    Nc_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    Ic_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    hw_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    bw_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    
    hb_b_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    hb_pho_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    hb_phi_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    hb_mix_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    coop_b_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    coop_pho_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)
    
    prot_counts=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype=np.int32)

    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            # build auxiliary boolean variable aa_present for given snapshot
            aa_belongs_to=np.zeros((box_size+1,box_size+1), dtype=np.int32)
            for i_prot in range(n_proteins):
                for i_mon in range(n_monomers):
                    aa_present[
                        all_x_pos[i_mon,i_prot,i_t,i_replica],
                        all_y_pos[i_mon,i_prot,i_t,i_replica],i_t,i_replica]=True
                    aa_belongs_to[
                        all_x_pos[i_mon,i_prot,i_t,i_replica],
                        all_y_pos[i_mon,i_prot,i_t,i_replica]]=i_prot
                    now_center_x=all_x_pos[i_mon,i_prot,i_t,i_replica]
                    now_center_y=all_y_pos[i_mon,i_prot,i_t,i_replica]
                    for i_x,i_y in [[-1,0],[1,0],[0,-1],[0,1]]:
                        # for i_y in [-1,1]:
                        now_x_ind=((now_center_x+i_x+box_size-1)%box_size)+1
                        now_y_ind=((now_center_y+i_y+box_size-1)%box_size)+1
                        hw_present[now_x_ind,now_y_ind,i_t,i_replica]=True
            for i_prot in range(n_proteins):
                for i_mon in range(n_monomers):
                    hw_present[all_x_pos[i_mon,i_prot,i_t,i_replica],
                        all_y_pos[i_mon,i_prot,i_t,i_replica],i_t,i_replica]=False

            # build auxiliary boolean variable aa_present for given snapshot
            for now_center_x in range(1,box_size+1):
                for now_center_y in range(1,box_size+1):
                    if (not(hw_present[now_center_x,now_center_y,i_t,i_replica]) and not(aa_present[now_center_x,now_center_y,i_t,i_replica])):
                        bw_present[now_center_x,now_center_y,i_t,i_replica]=True
                            
            # build auxiliary boolean variable aa_present for given snapshot
            for i_prot in range(n_proteins):
                for i_mon_A in range(n_monomers):
                    if nativ_neigh_list[i_mon_A]>=0:
                        i_mon_B=nativ_neigh_list[i_mon_A]
                        dx = box_size//2 - abs(abs(all_x_pos[i_mon_A,i_prot,i_t,i_replica]-all_x_pos[i_mon_B,i_prot,i_t,i_replica]) - box_size//2)
                        dy = box_size//2 - abs(abs(all_y_pos[i_mon_A,i_prot,i_t,i_replica]-all_y_pos[i_mon_B,i_prot,i_t,i_replica]) - box_size//2)
                        if ((dx == 1) and (dy == 0)) or ((dx == 0) and (dy == 1)):
                            Nc_present[all_x_pos[i_mon_A,i_prot,i_t,i_replica],all_y_pos[i_mon_A,i_prot,i_t,i_replica],i_t,i_replica]+=1
                            Nc_present[all_x_pos[i_mon_B,i_prot,i_t,i_replica],all_y_pos[i_mon_B,i_prot,i_t,i_replica],i_t,i_replica]+=1
           
            # build auxiliary boolean variable aa_present for given snapshot
            for i_prot_A in range(n_proteins-1):
                for i_prot_B in range( i_prot_A+1 , n_proteins):
                    for i_mon_A in range(n_monomers):
                        for i_mon_B in range(n_monomers):
                            dx = box_size//2 - abs(abs(all_x_pos[i_mon_A,i_prot_A,i_t,i_replica]-all_x_pos[i_mon_B,i_prot_B,i_t,i_replica]) - box_size//2)
                            dy = box_size//2 - abs(abs(all_y_pos[i_mon_A,i_prot_A,i_t,i_replica]-all_y_pos[i_mon_B,i_prot_B,i_t,i_replica]) - box_size//2)
                            if ((dx == 1) and (dy == 0)) or ((dx == 0) and (dy == 1)):
                                Ic_present[all_x_pos[i_mon_A,i_prot_A,i_t,i_replica],all_y_pos[i_mon_A,i_prot_A,i_t,i_replica],i_t,i_replica]+=1
                                Ic_present[all_x_pos[i_mon_B,i_prot_B,i_t,i_replica],all_y_pos[i_mon_B,i_prot_B,i_t,i_replica],i_t,i_replica]+=1
            
            # repeat counting for given snapshot for n_tests_per_snapshot times
            for now_center_x in range(1,box_size+1):
                for now_center_y in range(1,box_size+1):
                    now_prot_counts=np.zeros((n_cells_tested), dtype=np.int32)
                    now_prot_counts_ind=0
                    for i_x in range(-test_radius,test_radius+1):
                        for i_y in range(-test_radius,test_radius+1):
                            now_x_ind=((now_center_x+i_x+box_size-1)%box_size)+1
                            now_y_ind=((now_center_y+i_y+box_size-1)%box_size)+1
                            dx = box_size//2 - abs(abs(now_x_ind-now_center_x) - box_size//2)
                            dy = box_size//2 - abs(abs(now_y_ind-now_center_y) - box_size//2)
                            if (dx**2+dy**2 <= test_radius**2):
                                if aa_present[now_x_ind,now_y_ind,i_t,i_replica]:
                                    aa_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=1
                                    now_prot_counts[now_prot_counts_ind]=aa_belongs_to[now_x_ind,now_y_ind]
                                    now_prot_counts_ind+=1
                                if Nc_present[now_x_ind,now_y_ind,i_t,i_replica]>0:
                                    Nc_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=Nc_present[now_x_ind,now_y_ind,i_t,i_replica]
                                if Ic_present[now_x_ind,now_y_ind,i_t,i_replica]>0:
                                    Ic_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=Ic_present[now_x_ind,now_y_ind,i_t,i_replica]
                                if hw_present[now_x_ind,now_y_ind,i_t,i_replica]:
                                    hw_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=1
                                if bw_present[now_x_ind,now_y_ind,i_t,i_replica]:
                                    bw_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=1
                                if bulk_HB_present[now_x_ind,now_y_ind,i_t,i_replica]>0:
                                    hb_b_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=bulk_HB_present[now_x_ind,now_y_ind,i_t,i_replica]
                                if pho_HB_present[now_x_ind,now_y_ind,i_t,i_replica]>0:
                                    hb_pho_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=pho_HB_present[now_x_ind,now_y_ind,i_t,i_replica]
                                if phi_HB_present[now_x_ind,now_y_ind,i_t,i_replica]>0:
                                    hb_phi_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=phi_HB_present[now_x_ind,now_y_ind,i_t,i_replica]
                                if mix_HB_present[now_x_ind,now_y_ind,i_t,i_replica]>0:
                                    hb_mix_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=mix_HB_present[now_x_ind,now_y_ind,i_t,i_replica]
                                if coop_bonds_present[now_x_ind,now_y_ind,i_t,i_replica]>0:
                                    coop_b_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=coop_bonds_present[now_x_ind,now_y_ind,i_t,i_replica]
                                if coop_bonds_PHO_present[now_x_ind,now_y_ind,i_t,i_replica]>0:
                                    coop_pho_counts[now_center_x-1,now_center_y-1,i_t,i_replica]+=coop_bonds_PHO_present[now_x_ind,now_y_ind,i_t,i_replica]
                    prot_counts[now_center_x-1,now_center_y-1,i_t,i_replica]=len(np.unique(now_prot_counts))-1

    aa_present=aa_present[1:,1:,:,:].astype(np.uint8)
    Nc_present=Nc_present[1:,1:,:,:].astype(np.uint8)
    Ic_present=Ic_present[1:,1:,:,:].astype(np.uint8)   
    hw_present=hw_present[1:,1:,:,:].astype(np.uint8)   
    bw_present=bw_present[1:,1:,:,:].astype(np.uint8)
    aa_counts=aa_counts.astype(np.uint16)
    Nc_counts=Nc_counts.astype(np.uint16)
    Ic_counts=Ic_counts.astype(np.uint16)
    hw_counts=hw_counts.astype(np.uint16)
    bw_counts=bw_counts.astype(np.uint16)
    hb_b_counts=hb_b_counts.astype(np.uint16)
    hb_pho_counts=hb_pho_counts.astype(np.uint16)
    hb_phi_counts=hb_phi_counts.astype(np.uint16)
    hb_mix_counts=hb_mix_counts.astype(np.uint16)
    coop_b_counts=coop_b_counts.astype(np.uint16)
    coop_pho_counts=coop_pho_counts.astype(np.uint16)
    
    prot_counts=prot_counts.astype(np.uint16)

    bulk_HB_present=bulk_HB_present[1:,1:,:,:].astype(np.uint8)
    pho_HB_present=pho_HB_present[1:,1:,:,:].astype(np.uint8)
    phi_HB_present=phi_HB_present[1:,1:,:,:].astype(np.uint8)
    mix_HB_present=mix_HB_present[1:,1:,:,:].astype(np.uint8)
    coop_bonds_present=coop_bonds_present[1:,1:,:,:].astype(np.uint8)
    coop_bonds_PHO_present=coop_bonds_PHO_present[1:,1:,:,:].astype(np.uint8)
    
    
    return aa_present,Nc_present,Ic_present,hw_present,bw_present,\
        aa_counts,Nc_counts,Ic_counts,hw_counts,bw_counts,\
        hb_b_counts,hb_pho_counts,hb_phi_counts,hb_mix_counts,coop_b_counts,coop_pho_counts,\
        bulk_HB_present,pho_HB_present,phi_HB_present,mix_HB_present,coop_bonds_present,coop_bonds_PHO_present,\
        prot_counts,\
        n_cells_tested





@njit
def build_aux_vars_hist_2d(local_counts_Nc,n_proteins,n_native_bonds_per_protein,n_cells_tested):
    aux_var_1=np.zeros(local_counts_Nc.sum(),dtype=np.int32)
    aux_var_2=np.zeros(local_counts_Nc.sum(),dtype=np.int32)
    ii=0
    for iNc in range(n_proteins*n_native_bonds_per_protein+1):
        for iCounts in range(n_cells_tested+1):
            for _ in range(local_counts_Nc[iNc,iCounts]):
                aux_var_1[ii]=iNc
                aux_var_2[ii]=iCounts
                ii+=1
    return aux_var_1,aux_var_2


def retrieve_simul_variables(all_x_pos,nativ_neigh_list):
    n_monomers=all_x_pos.shape[0]
    n_proteins=all_x_pos.shape[1]
    n_snapshots=all_x_pos.shape[2]
    n_replicas=all_x_pos.shape[3]
    n_native_bonds_per_protein=np.sum(nativ_neigh_list>0)
    return n_monomers,n_proteins,n_snapshots,n_replicas,n_native_bonds_per_protein



    
def catenate_files_local(npz_loaded_file_list):
    if type(npz_loaded_file_list) == list:
        box_size=int(npz_loaded_file_list[0]['box_size'])
        nativ_neigh_list=npz_loaded_file_list[0]['nativ_neigh_list']
        path_to_array_replica_folder=npz_loaded_file_list[0]['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list[0]['all_x_pos']
        all_y_pos=npz_loaded_file_list[0]['all_y_pos']
        all_Nc=npz_loaded_file_list[0]['all_Nc']
        all_Ic=npz_loaded_file_list[0]['all_Ic']
        list_of_times=npz_loaded_file_list[0]['list_of_times']
        all_Nc_per_prot=npz_loaded_file_list[0]['all_Nc_per_prot']
        aa_counts=npz_loaded_file_list[0]['aa_counts']
        Nc_counts=npz_loaded_file_list[0]['Nc_counts']
        Ic_counts=npz_loaded_file_list[0]['Ic_counts']
        hw_counts=npz_loaded_file_list[0]['hw_counts']
        bw_counts=npz_loaded_file_list[0]['bw_counts']
        n_cells_tested=npz_loaded_file_list[0]['n_cells_tested']
        raw_data_folding=npz_loaded_file_list[0]['raw_data_folding']            
        
        for i_file in range(1,len(npz_loaded_file_list)):
            all_x_pos=np.append(all_x_pos,npz_loaded_file_list[i_file]['all_x_pos'],axis=2)
            all_y_pos=np.append(all_y_pos,npz_loaded_file_list[i_file]['all_y_pos'],axis=2)
            all_Nc=np.append(all_Nc,npz_loaded_file_list[i_file]['all_Nc'],axis=1)
            all_Ic=np.append(all_Ic,npz_loaded_file_list[i_file]['all_Ic'],axis=1)
            list_of_times=np.append(list_of_times,npz_loaded_file_list[i_file]['list_of_times']+list_of_times[-1])
            all_Nc_per_prot=np.append(all_Nc_per_prot,npz_loaded_file_list[i_file]['all_Nc_per_prot'],axis=2)
            aa_counts=np.append(aa_counts,npz_loaded_file_list[i_file]['aa_counts'],axis=2)
            Nc_counts=np.append(Nc_counts,npz_loaded_file_list[i_file]['Nc_counts'],axis=2)
            Ic_counts=np.append(Ic_counts,npz_loaded_file_list[i_file]['Ic_counts'],axis=2)
            hw_counts=np.append(hw_counts,npz_loaded_file_list[i_file]['hw_counts'],axis=2)
            bw_counts=np.append(bw_counts,npz_loaded_file_list[i_file]['bw_counts'],axis=2)
            raw_data_folding=np.append(raw_data_folding,npz_loaded_file_list[i_file]['raw_data_folding'],axis=0)
    else:
        box_size=int(npz_loaded_file_list['box_size'])
        nativ_neigh_list=npz_loaded_file_list['nativ_neigh_list']
        path_to_array_replica_folder=npz_loaded_file_list['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list['all_x_pos']
        all_y_pos=npz_loaded_file_list['all_y_pos']
        all_Nc=npz_loaded_file_list['all_Nc']
        all_Ic=npz_loaded_file_list['all_Ic']
        list_of_times=npz_loaded_file_list['list_of_times']
        all_Nc_per_prot=npz_loaded_file_list['all_Nc_per_prot']
        aa_counts=npz_loaded_file_list['aa_counts']
        Nc_counts=npz_loaded_file_list['Nc_counts']
        Ic_counts=npz_loaded_file_list['Ic_counts']
        hw_counts=npz_loaded_file_list['hw_counts']
        bw_counts=npz_loaded_file_list['bw_counts']
        n_cells_tested=npz_loaded_file_list['n_cells_tested']
        raw_data_folding=npz_loaded_file_list['raw_data_folding']
        
    return box_size,nativ_neigh_list,path_to_array_replica_folder,\
    all_x_pos,all_y_pos,all_Nc,all_Ic,list_of_times,all_Nc_per_prot,\
    aa_counts,Nc_counts,Ic_counts,hw_counts,bw_counts,n_cells_tested,\
    raw_data_folding



    
def catenate_files_all_thermo(npz_loaded_file_list):
    if type(npz_loaded_file_list) == list:
        box_size=int(npz_loaded_file_list[0]['box_size'])
        nativ_neigh_list=npz_loaded_file_list[0]['nativ_neigh_list']
        path_to_array_replica_folder=npz_loaded_file_list[0]['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list[0]['all_x_pos']
        all_y_pos=npz_loaded_file_list[0]['all_y_pos']
        all_Nc=npz_loaded_file_list[0]['all_Nc']
        all_Ic=npz_loaded_file_list[0]['all_Ic']
        list_of_times=npz_loaded_file_list[0]['list_of_times']
        all_Nc_per_prot=npz_loaded_file_list[0]['all_Nc_per_prot']
        bulk_molecules_rad_10=npz_loaded_file_list[0]['bulk_molecules_rad_10']
        water_molecules_10=npz_loaded_file_list[0]['water_molecules_10']
        bulk_molecules_rad_5=npz_loaded_file_list[0]['bulk_molecules_rad_5']
        water_molecules_5=npz_loaded_file_list[0]['water_molecules_5']
        local_counts_frequencies_5=npz_loaded_file_list[0]['local_counts_frequencies_5']
        local_counts_Nc_5=npz_loaded_file_list[0]['local_counts_Nc_5']
        local_counts_Ic_5=npz_loaded_file_list[0]['local_counts_Ic_5']
        n_cells_tested_5=npz_loaded_file_list[0]['n_cells_tested_5']
        local_counts_frequencies_10=npz_loaded_file_list[0]['local_counts_frequencies_10']
        local_counts_Nc_10=npz_loaded_file_list[0]['local_counts_Nc_10']
        local_counts_Ic_10=npz_loaded_file_list[0]['local_counts_Ic_10']
        n_cells_tested_10=npz_loaded_file_list[0]['n_cells_tested_10']
        raw_data_folding=npz_loaded_file_list[0]['raw_data_folding']
        for i_file in range(1,len(npz_loaded_file_list)):
            all_x_pos=np.append(all_x_pos,npz_loaded_file_list[i_file]['all_x_pos'],axis=2)
            all_y_pos=np.append(all_y_pos,npz_loaded_file_list[i_file]['all_y_pos'],axis=2)
            all_Nc=np.append(all_Nc,npz_loaded_file_list[i_file]['all_Nc'],axis=1)
            all_Ic=np.append(all_Ic,npz_loaded_file_list[i_file]['all_Ic'],axis=1)
            list_of_times=np.append(list_of_times,npz_loaded_file_list[i_file]['list_of_times']+list_of_times[-1])
            all_Nc_per_prot=np.append(all_Nc_per_prot,npz_loaded_file_list[i_file]['all_Nc_per_prot'],axis=2)
            bulk_molecules_rad_10=np.append(bulk_molecules_rad_10,npz_loaded_file_list[i_file]['bulk_molecules_rad_10'],axis=1)
            bulk_molecules_rad_5=np.append(bulk_molecules_rad_5,npz_loaded_file_list[i_file]['bulk_molecules_rad_5'],axis=1)
            local_counts_frequencies_5=np.append(local_counts_frequencies_5,npz_loaded_file_list[i_file]['local_counts_frequencies_5'],axis=1)
            local_counts_frequencies_10=np.append(local_counts_frequencies_10,npz_loaded_file_list[i_file]['local_counts_frequencies_10'],axis=1)
            local_counts_Nc_5=local_counts_Nc_5+npz_loaded_file_list[i_file]['local_counts_Nc_5']
            local_counts_Ic_5=local_counts_Ic_5+npz_loaded_file_list[i_file]['local_counts_Ic_5']
            local_counts_Nc_10=local_counts_Nc_10+npz_loaded_file_list[i_file]['local_counts_Nc_10']
            local_counts_Ic_10=local_counts_Ic_10+npz_loaded_file_list[i_file]['local_counts_Ic_10']
            raw_data_folding=np.append(raw_data_folding,npz_loaded_file_list[i_file]['raw_data_folding'],axis=0)
    else:
        box_size=int(npz_loaded_file_list['box_size'])
        nativ_neigh_list=npz_loaded_file_list['nativ_neigh_list']
        path_to_array_replica_folder=npz_loaded_file_list['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list['all_x_pos']
        all_y_pos=npz_loaded_file_list['all_y_pos']
        all_Nc=npz_loaded_file_list['all_Nc']
        all_Ic=npz_loaded_file_list['all_Ic']
        list_of_times=npz_loaded_file_list['list_of_times']
        all_Nc_per_prot=npz_loaded_file_list['all_Nc_per_prot']
        bulk_molecules_rad_10=npz_loaded_file_list['bulk_molecules_rad_10']
        water_molecules_10=npz_loaded_file_list['water_molecules_10']
        bulk_molecules_rad_5=npz_loaded_file_list['bulk_molecules_rad_5']
        water_molecules_5=npz_loaded_file_list['water_molecules_5']
        local_counts_frequencies_5=npz_loaded_file_list['local_counts_frequencies_5']
        local_counts_Nc_5=npz_loaded_file_list['local_counts_Nc_5']
        local_counts_Ic_5=npz_loaded_file_list['local_counts_Ic_5']
        n_cells_tested_5=npz_loaded_file_list['n_cells_tested_5']
        local_counts_frequencies_10=npz_loaded_file_list['local_counts_frequencies_10']
        local_counts_Nc_10=npz_loaded_file_list['local_counts_Nc_10']
        local_counts_Ic_10=npz_loaded_file_list['local_counts_Ic_10']
        n_cells_tested_10=npz_loaded_file_list['n_cells_tested_10']
        raw_data_folding=npz_loaded_file_list['raw_data_folding']
        
    return box_size,nativ_neigh_list,path_to_array_replica_folder,\
    all_x_pos,all_y_pos,all_Nc,all_Ic,list_of_times,\
    all_Nc_per_prot,\
    bulk_molecules_rad_10,water_molecules_10,bulk_molecules_rad_5,water_molecules_5,\
    local_counts_frequencies_5,local_counts_Nc_5,local_counts_Ic_5,n_cells_tested_5,\
    local_counts_frequencies_10,local_counts_Nc_10,local_counts_Ic_10,n_cells_tested_10,\
    raw_data_folding

    
def catenate_files_all(npz_loaded_file_list):
    if type(npz_loaded_file_list) == list:
        box_size=int(npz_loaded_file_list[0]['box_size'])
        nativ_neigh_list=npz_loaded_file_list[0]['nativ_neigh_list']
        path_to_array_replica_folder=npz_loaded_file_list[0]['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list[0]['all_x_pos']
        all_y_pos=npz_loaded_file_list[0]['all_y_pos']
        all_Nc=npz_loaded_file_list[0]['all_Nc']
        all_Ic=npz_loaded_file_list[0]['all_Ic']
        list_of_times=npz_loaded_file_list[0]['list_of_times']
        all_Nc_per_prot=npz_loaded_file_list[0]['all_Nc_per_prot']
        bulk_molecules_rad_10=npz_loaded_file_list[0]['bulk_molecules_rad_10']
        water_molecules_10=npz_loaded_file_list[0]['water_molecules_10']
        bulk_molecules_rad_5=npz_loaded_file_list[0]['bulk_molecules_rad_5']
        water_molecules_5=npz_loaded_file_list[0]['water_molecules_5']
        local_counts_frequencies_5=npz_loaded_file_list[0]['local_counts_frequencies_5']
        local_counts_Nc_5=npz_loaded_file_list[0]['local_counts_Nc_5']
        local_counts_Ic_5=npz_loaded_file_list[0]['local_counts_Ic_5']
        n_cells_tested_5=npz_loaded_file_list[0]['n_cells_tested_5']
        local_counts_frequencies_10=npz_loaded_file_list[0]['local_counts_frequencies_10']
        local_counts_Nc_10=npz_loaded_file_list[0]['local_counts_Nc_10']
        local_counts_Ic_10=npz_loaded_file_list[0]['local_counts_Ic_10']
        n_cells_tested_10=npz_loaded_file_list[0]['n_cells_tested_10']
        for i_file in range(1,len(npz_loaded_file_list)):
            all_x_pos=np.append(all_x_pos,npz_loaded_file_list[i_file]['all_x_pos'],axis=2)
            all_y_pos=np.append(all_y_pos,npz_loaded_file_list[i_file]['all_y_pos'],axis=2)
            all_Nc=np.append(all_Nc,npz_loaded_file_list[i_file]['all_Nc'],axis=1)
            all_Ic=np.append(all_Ic,npz_loaded_file_list[i_file]['all_Ic'],axis=1)
            list_of_times=np.append(list_of_times,npz_loaded_file_list[i_file]['list_of_times']+list_of_times[-1])
            all_Nc_per_prot=np.append(all_Nc_per_prot,npz_loaded_file_list[i_file]['all_Nc_per_prot'],axis=2)
            bulk_molecules_rad_10=np.append(bulk_molecules_rad_10,npz_loaded_file_list[i_file]['bulk_molecules_rad_10'],axis=1)
            bulk_molecules_rad_5=np.append(bulk_molecules_rad_5,npz_loaded_file_list[i_file]['bulk_molecules_rad_5'],axis=1)
            local_counts_frequencies_5=np.append(local_counts_frequencies_5,npz_loaded_file_list[i_file]['local_counts_frequencies_5'],axis=1)
            local_counts_frequencies_10=np.append(local_counts_frequencies_10,npz_loaded_file_list[i_file]['local_counts_frequencies_10'],axis=1)
            local_counts_Nc_5=local_counts_Nc_5+npz_loaded_file_list[i_file]['local_counts_Nc_5']
            local_counts_Ic_5=local_counts_Ic_5+npz_loaded_file_list[i_file]['local_counts_Ic_5']
            local_counts_Nc_10=local_counts_Nc_10+npz_loaded_file_list[i_file]['local_counts_Nc_10']
            local_counts_Ic_10=local_counts_Ic_10+npz_loaded_file_list[i_file]['local_counts_Ic_10']
    else:
        box_size=int(npz_loaded_file_list['box_size'])
        nativ_neigh_list=npz_loaded_file_list['nativ_neigh_list']
        path_to_array_replica_folder=npz_loaded_file_list['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list['all_x_pos']
        all_y_pos=npz_loaded_file_list['all_y_pos']
        all_Nc=npz_loaded_file_list['all_Nc']
        all_Ic=npz_loaded_file_list['all_Ic']
        list_of_times=npz_loaded_file_list['list_of_times']
        all_Nc_per_prot=npz_loaded_file_list['all_Nc_per_prot']
        bulk_molecules_rad_10=npz_loaded_file_list['bulk_molecules_rad_10']
        water_molecules_10=npz_loaded_file_list['water_molecules_10']
        bulk_molecules_rad_5=npz_loaded_file_list['bulk_molecules_rad_5']
        water_molecules_5=npz_loaded_file_list['water_molecules_5']
        local_counts_frequencies_5=npz_loaded_file_list['local_counts_frequencies_5']
        local_counts_Nc_5=npz_loaded_file_list['local_counts_Nc_5']
        local_counts_Ic_5=npz_loaded_file_list['local_counts_Ic_5']
        n_cells_tested_5=npz_loaded_file_list['n_cells_tested_5']
        local_counts_frequencies_10=npz_loaded_file_list['local_counts_frequencies_10']
        local_counts_Nc_10=npz_loaded_file_list['local_counts_Nc_10']
        local_counts_Ic_10=npz_loaded_file_list['local_counts_Ic_10']
        n_cells_tested_10=npz_loaded_file_list['n_cells_tested_10']
        
    return box_size,nativ_neigh_list,path_to_array_replica_folder,\
    all_x_pos,all_y_pos,all_Nc,all_Ic,list_of_times,\
    all_Nc_per_prot,\
    bulk_molecules_rad_10,water_molecules_10,bulk_molecules_rad_5,water_molecules_5,\
    local_counts_frequencies_5,local_counts_Nc_5,local_counts_Ic_5,n_cells_tested_5,\
    local_counts_frequencies_10,local_counts_Nc_10,local_counts_Ic_10,n_cells_tested_10



@njit
def calculate_all_local_concentration_systematic(all_x_pos,all_y_pos,box_size,test_radius,all_Nc,all_Ic,n_native_bonds_per_protein):
    n_monomers=all_x_pos.shape[0]
    n_proteins=all_x_pos.shape[1]
    n_snapshots=all_x_pos.shape[2]
    n_replicas=all_x_pos.shape[3]
    
    # count number of cells involved in calcualtion for given test_radius
    n_cells_tested=0
    now_center_x=1
    now_center_y=1
    for i_x in range(-test_radius,test_radius+1):
        for i_y in range(-test_radius,test_radius+1):
            now_x_ind=((now_center_x+i_x+box_size-1)%box_size)+1
            now_y_ind=((now_center_y+i_y+box_size-1)%box_size)+1
            dx = box_size//2 - abs(abs(now_x_ind-now_center_x) - box_size//2)
            dy = box_size//2 - abs(abs(now_y_ind-now_center_y) - box_size//2)
            if (dx**2+dy**2 <= test_radius**2):
                n_cells_tested+=1
    # count how many cells inside radius are filled
    # repeat for every replica and every snapshot
    #
    local_counts_frequencies=np.zeros((n_replicas,n_snapshots,n_cells_tested+1), dtype=np.int32)
    local_counts_Nc=np.zeros((n_proteins*n_native_bonds_per_protein+1,n_cells_tested+1), dtype=np.int32)
    local_counts_Ic=np.zeros((n_proteins*n_monomers+1,n_cells_tested+1), dtype=np.int32)
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            # build auxiliary boolean variable system_state for given snapshot
            system_state=np.zeros((box_size+1,box_size+1), dtype=np.bool_)
            for i_prot in range(n_proteins):
                for i_mon in range(n_monomers):
                    system_state[
                        all_x_pos[i_mon,i_prot,i_t,i_replica],
                        all_y_pos[i_mon,i_prot,i_t,i_replica]]=True
            # repeat counting for given snapshot for n_tests_per_snapshot times
            for now_center_x in range(1,box_size+1):
                for now_center_y in range(1,box_size+1):
                    now_counts=0
                    for i_x in range(-test_radius,test_radius+1):
                        for i_y in range(-test_radius,test_radius+1):
                            now_x_ind=((now_center_x+i_x+box_size-1)%box_size)+1
                            now_y_ind=((now_center_y+i_y+box_size-1)%box_size)+1
                            dx = box_size//2 - abs(abs(now_x_ind-now_center_x) - box_size//2)
                            dy = box_size//2 - abs(abs(now_y_ind-now_center_y) - box_size//2)
                            if ((dx**2+dy**2 <= test_radius**2) and system_state[now_x_ind,now_y_ind]) :
                                now_counts+=1
                    local_counts_frequencies[i_replica,i_t,now_counts]+=1
                    local_counts_Nc[all_Nc[i_replica,i_t],now_counts]+=1
                    local_counts_Ic[all_Ic[i_replica,i_t],now_counts]+=1
    
    return local_counts_frequencies,local_counts_Nc,local_counts_Ic,n_cells_tested






@njit
def count_fold_unfold_prot(all_x_pos,all_y_pos,box_size,nativ_neigh_list):
    n_monomers=all_x_pos.shape[0]
    n_proteins=all_x_pos.shape[1]
    n_snapshots=all_x_pos.shape[2]
    n_replicas=all_x_pos.shape[3]
    
    all_Nc_per_prot=np.zeros((n_proteins,n_replicas,n_snapshots),dtype=np.int32)
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            for i_prot in range(n_proteins):
                for i_mon_A in range(n_monomers):
                    if nativ_neigh_list[i_mon_A]>=0:
                        i_mon_B=nativ_neigh_list[i_mon_A]
                        dx = box_size//2 - abs(abs(all_x_pos[i_mon_A,i_prot,i_t,i_replica]-all_x_pos[i_mon_B,i_prot,i_t,i_replica]) - box_size//2)
                        dy = box_size//2 - abs(abs(all_y_pos[i_mon_A,i_prot,i_t,i_replica]-all_y_pos[i_mon_B,i_prot,i_t,i_replica]) - box_size//2)
                        if ((dx == 1) and (dy == 0)) or ((dx == 0) and (dy == 1)):
                            all_Nc_per_prot[i_prot,i_replica,i_t]+=1
    return all_Nc_per_prot


@njit
def calculate_all_bulk_fractions_radius(all_x_pos,all_y_pos,box_size,test_radius):
    n_monomers=all_x_pos.shape[0]
    n_proteins=all_x_pos.shape[1]
    n_snapshots=all_x_pos.shape[2]
    n_replicas=all_x_pos.shape[3]    
    bulk_molecules=np.zeros((n_replicas,n_snapshots), dtype=np.int32)
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            system_state=np.zeros((box_size,box_size), dtype=np.bool_)
            for i_prot in range(n_proteins):
                for i_mon in range(n_monomers):
                    system_state[
                        all_x_pos[i_mon,i_prot,i_t,i_replica]-1,
                        all_y_pos[i_mon,i_prot,i_t,i_replica]-1]=True
            water_molecules=0
            for now_center_x in range(box_size):
                for now_center_y in range(box_size):
                    if not(system_state[now_center_x,now_center_y]):
                        water_molecules+=1
                        there_is_no_aa_near=True
                        for i_x in range(-test_radius,test_radius+1):
                            for i_y in range(-test_radius,test_radius+1):
                                now_x_ind=((now_center_x+i_x+box_size)%box_size)
                                now_y_ind=((now_center_y+i_y+box_size)%box_size)
                                dx = abs(now_x_ind-now_center_x)
                                if dx > box_size/2: dx=box_size-dx
                                dy = abs(now_y_ind-now_center_y)
                                if dy > box_size/2: dy=box_size-dy
                                if (dx**2+dy**2 <= test_radius**2):
                                    if system_state[now_x_ind,now_y_ind]:
                                        there_is_no_aa_near=False
                                if not(there_is_no_aa_near): break
                            if not(there_is_no_aa_near): break
                        if there_is_no_aa_near:
                            bulk_molecules[i_replica,i_t]+=1
    return bulk_molecules,water_molecules

@njit
def calculate_all_bulk_fractions(all_x_pos,all_y_pos,box_size):
    n_monomers=all_x_pos.shape[0]
    n_proteins=all_x_pos.shape[1]
    n_snapshots=all_x_pos.shape[2]
    n_replicas=all_x_pos.shape[3]    
    bulk_molecules=np.zeros((n_replicas,n_snapshots), dtype=np.int32)
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            system_state=np.zeros((box_size,box_size), dtype=np.bool_)
            for i_prot in range(n_proteins):
                for i_mon in range(n_monomers):
                    system_state[
                        all_x_pos[i_mon,i_prot,i_t,i_replica]-1,
                        all_y_pos[i_mon,i_prot,i_t,i_replica]-1]=True
            water_molecules=0
            for i_x in range(box_size):
                now_ln=((i_x+box_size-1)%box_size)
                now_rn=((i_x+box_size+1)%box_size)
                for i_y in range(box_size):
                    now_dn=((i_y+box_size-1)%box_size)
                    now_un=((i_y+box_size+1)%box_size)
                    if not(system_state[i_x,i_y]):
                        water_molecules+=1
                        if (not(system_state[i_x,now_dn]) and 
                            not(system_state[i_x,now_un]) and 
                            not(system_state[now_ln,i_y]) and 
                            not(system_state[now_rn,i_y])):
                            bulk_molecules[i_replica,i_t]+=1
    return bulk_molecules,water_molecules


@njit
def calculate_all_local_concentration(all_x_pos,all_y_pos,box_size,test_radius,n_tests_per_snapshot,all_Nc,all_Ic,n_native_bonds_per_protein):
    n_monomers=all_x_pos.shape[0]
    n_proteins=all_x_pos.shape[1]
    n_snapshots=all_x_pos.shape[2]
    n_replicas=all_x_pos.shape[3]
    
    # count number of cells involved in calcualtion for given test_radius
    n_cells_tested=0
    now_center_x=int((box_size*np.random.rand()+1)//1)
    now_center_y=int((box_size*np.random.rand()+1)//1)
    for i_x in range(-test_radius,test_radius+1):
        for i_y in range(-test_radius,test_radius+1):
            now_x_ind=((now_center_x+i_x+box_size-1)%box_size)+1
            now_y_ind=((now_center_y+i_y+box_size-1)%box_size)+1
            dx = box_size//2 - abs(abs(now_x_ind-now_center_x) - box_size//2)
            dy = box_size//2 - abs(abs(now_y_ind-now_center_y) - box_size//2)
            if (dx**2+dy**2 <= test_radius**2):
                n_cells_tested+=1
    
    # count how many cells inside radius are filled
    # repeat for every replica and every snapshot
    #
    local_counts_frequencies=np.zeros((n_replicas,n_snapshots,n_cells_tested), dtype=np.int32)
    local_counts_Nc=np.zeros((n_proteins*n_native_bonds_per_protein,n_cells_tested), dtype=np.int32)
    local_counts_Ic=np.zeros((n_proteins*n_monomers,n_cells_tested), dtype=np.int32)
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            # build auxiliary boolean variable system_state for given snapshot
            system_state=np.zeros((box_size+1,box_size+1), dtype=np.bool_)
            for i_prot in range(n_proteins):
                for i_mon in range(n_monomers):
                    system_state[
                        all_x_pos[i_mon,i_prot,i_t,i_replica],
                        all_y_pos[i_mon,i_prot,i_t,i_replica]]=True
            # repeat counting for given snapshot for n_tests_per_snapshot times
            for i_test in range(n_tests_per_snapshot):
                # extract randomly circle center:
                now_center_x=int((box_size*np.random.rand()+1)//1)
                now_center_y=int((box_size*np.random.rand()+1)//1)
                # extract randomly circle center:
                now_counts=0
                for i_x in range(-test_radius,test_radius+1):
                    for i_y in range(-test_radius,test_radius+1):
                        now_x_ind=((now_center_x+i_x+box_size-1)%box_size)+1
                        now_y_ind=((now_center_y+i_y+box_size-1)%box_size)+1
                        dx = box_size//2 - abs(abs(now_x_ind-now_center_x) - box_size//2)
                        dy = box_size//2 - abs(abs(now_y_ind-now_center_y) - box_size//2)
                        if ((dx**2+dy**2 <= test_radius**2) and system_state[now_x_ind,now_y_ind]) :
                            now_counts+=1
                local_counts_frequencies[i_replica,i_t,now_counts]+=1
                local_counts_Nc[all_Nc[i_replica,i_t],now_counts]+=1
                local_counts_Ic[all_Ic[i_replica,i_t],now_counts]+=1
    
    return local_counts_frequencies,local_counts_Nc,local_counts_Ic,n_cells_tested

def reconstruct_variables(npz_loaded_file):
    box_size=int(npz_loaded_file['L'])
    n_proteins=npz_loaded_file['n_proteins']
    n_replicas=npz_loaded_file['n_replicas']
    nativ_neigh_list=npz_loaded_file['nativ_neigh_list']
    n_native_bonds_per_protein=npz_loaded_file['n_native_bonds_per_protein']
    path_to_array_replica_folder=npz_loaded_file['path_to_array_replica_folder']
    all_x_pos=npz_loaded_file['all_x_pos']
    all_y_pos=npz_loaded_file['all_y_pos']
    all_Nc=npz_loaded_file['all_Nc']
    all_Ic=npz_loaded_file['all_Ic']
    list_of_times=npz_loaded_file['list_of_times']
    return box_size,n_proteins,n_replicas,nativ_neigh_list,n_native_bonds_per_protein,path_to_array_replica_folder, \
    all_x_pos,all_y_pos,all_Nc,all_Ic,list_of_times
    
def catenate_files(npz_loaded_file_list):
    if type(npz_loaded_file_list) == list:
        box_size=int(npz_loaded_file_list[0]['L'])
        n_proteins=npz_loaded_file_list[0]['n_proteins']
        n_replicas=npz_loaded_file_list[0]['n_replicas']
        nativ_neigh_list=npz_loaded_file_list[0]['nativ_neigh_list']
        n_native_bonds_per_protein=npz_loaded_file_list[0]['n_native_bonds_per_protein']
        path_to_array_replica_folder=npz_loaded_file_list[0]['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list[0]['all_x_pos']
        all_y_pos=npz_loaded_file_list[0]['all_y_pos']
        all_Nc=npz_loaded_file_list[0]['all_Nc']
        all_Ic=npz_loaded_file_list[0]['all_Ic']
        list_of_times=npz_loaded_file_list[0]['list_of_times']
        for i_file in range(1,len(npz_loaded_file_list)):
            all_x_pos=np.append(all_x_pos,npz_loaded_file_list[i_file]['all_x_pos'],axis=2)
            all_y_pos=np.append(all_y_pos,npz_loaded_file_list[i_file]['all_y_pos'],axis=2)
            all_Nc=np.append(all_Nc,npz_loaded_file_list[i_file]['all_Nc'],axis=1)
            all_Ic=np.append(all_Ic,npz_loaded_file_list[i_file]['all_Ic'],axis=1)
            list_of_times=np.append(list_of_times,npz_loaded_file_list[i_file]['list_of_times']+list_of_times[-1])
    else:
        box_size=int(npz_loaded_file_list['L'])
        n_proteins=npz_loaded_file_list['n_proteins']
        n_replicas=npz_loaded_file_list['n_replicas']
        nativ_neigh_list=npz_loaded_file_list['nativ_neigh_list']
        n_native_bonds_per_protein=npz_loaded_file_list['n_native_bonds_per_protein']
        path_to_array_replica_folder=npz_loaded_file_list['path_to_array_replica_folder']
        all_x_pos=npz_loaded_file_list['all_x_pos']
        all_y_pos=npz_loaded_file_list['all_y_pos']
        all_Nc=npz_loaded_file_list['all_Nc']
        all_Ic=npz_loaded_file_list['all_Ic']
        list_of_times=npz_loaded_file_list['list_of_times']
        
    return box_size,n_proteins,n_replicas,nativ_neigh_list,n_native_bonds_per_protein,path_to_array_replica_folder, \
    all_x_pos,all_y_pos,all_Nc,all_Ic,list_of_times

def read_conf_file(path,n_proteins,n_monomers,max_n_snapshots):
    raw_conf_data=np.loadtxt(path,dtype=np.int32,max_rows=max_n_snapshots*n_proteins*n_monomers)
    return raw_conf_data


def convert_raw_conf_data(raw_conf_data,n_proteins,n_monomers):
    raw_t_data=raw_conf_data[:,0]
    raw_x_data=raw_conf_data[:,1]
    raw_y_data=raw_conf_data[:,2]

    list_of_times=np.unique(raw_t_data)
    t_i=list_of_times[0]
    t_f=list_of_times[-1]
    d_t=list_of_times[1]-list_of_times[0]

    n_snapshots=list_of_times.size
    x_pos=np.empty([n_monomers, n_proteins, n_snapshots],dtype=np.int32)
    y_pos=np.empty([n_monomers, n_proteins, n_snapshots],dtype=np.int32)

    for i_t in range(n_snapshots):
        now_t_indexes=(raw_t_data==list_of_times[i_t])
        now_x=raw_x_data[now_t_indexes]
        now_y=raw_y_data[now_t_indexes]
        for i_prot in range(n_proteins):
            now_prot_indexes = np.arange(0,n_monomers) + i_prot*n_monomers
            x_pos[:,i_prot,i_t]=now_x[now_prot_indexes]
            y_pos[:,i_prot,i_t]=now_y[now_prot_indexes]
    return x_pos,y_pos,list_of_times

def read_target_structures(path,n_monomers):
    native_structure=np.loadtxt(path,dtype=np.int32,max_rows=n_monomers)
    nativ_neigh_list=np.ones(n_monomers,dtype=np.int32)*-1
    for i_mon_A in range(0,n_monomers):
        for i_mon_B in range(i_mon_A+2,n_monomers):
            dx = abs(native_structure[i_mon_A,0]-native_structure[i_mon_B,0])
            dy = abs(native_structure[i_mon_A,1]-native_structure[i_mon_B,1])
            if ((dx == 1) and (dy == 0)) or ((dx == 0) and (dy == 1)):
                nativ_neigh_list[i_mon_A]=i_mon_B
    return nativ_neigh_list

@njit
def count_all_native_neigh(x_pos,y_pos,nativ_neigh_list,box_size):
    n_monomers=x_pos.shape[0]
    n_proteins=x_pos.shape[1]
    n_snapshots=x_pos.shape[2]
    all_Nc=np.zeros(n_snapshots,dtype=np.int32)
    for i_t in range(n_snapshots):
        for i_prot in range(n_proteins):
            for i_mon_A in range(n_monomers):
                if nativ_neigh_list[i_mon_A]>=0:
                    i_mon_B=nativ_neigh_list[i_mon_A]
                    dx = box_size//2 - abs(abs(x_pos[i_mon_A,i_prot,i_t]-x_pos[i_mon_B,i_prot,i_t]) - box_size//2)
                    dy = box_size//2 - abs(abs(y_pos[i_mon_A,i_prot,i_t]-y_pos[i_mon_B,i_prot,i_t]) - box_size//2)
                    if ((dx == 1) and (dy == 0)) or ((dx == 0) and (dy == 1)):
                        all_Nc[i_t]+=1
    return all_Nc

@njit
def count_all_inter_contacts(x_pos,y_pos,box_size):
    n_monomers=x_pos.shape[0]
    n_proteins=x_pos.shape[1]
    n_snapshots=x_pos.shape[2]
    all_Ic=np.zeros(n_snapshots, dtype=np.int32)
    for i_t in range(n_snapshots):
        for i_prot_A in range(n_proteins-1):
            for i_prot_B in range( i_prot_A+1 , n_proteins):
                for i_mon_A in range(n_monomers):
                    for i_mon_B in range(n_monomers):
                        dx = box_size//2 - abs(abs(x_pos[i_mon_A,i_prot_A,i_t]-x_pos[i_mon_B,i_prot_B,i_t]) - box_size//2)
                        dy = box_size//2 - abs(abs(y_pos[i_mon_A,i_prot_A,i_t]-y_pos[i_mon_B,i_prot_B,i_t]) - box_size//2)
                        if ((dx == 1) and (dy == 0)) or ((dx == 0) and (dy == 1)):
                            all_Ic[i_t] += 1
    return all_Ic

@njit
def correlation_function(time_series):
    time_series_mean=np.mean(time_series)
    T=len(time_series)
    Ck=np.zeros(T//2)
    for lag in range(T//2):
        cnt = 0
        accum = 0.0
        for i in range(T-lag):
            accum += (time_series[i]-time_series_mean)*(time_series[i+lag]-time_series_mean)
            cnt += 1
        Ck[lag] = accum/cnt
    return Ck







# %% definitions of cluster functions
@njit
def calculate_radius_of_gyration_clusters(liq_cluster_label, liq_cluster_size):
    box_size=liq_cluster_label.shape[0]
    n_snapshots=liq_cluster_label.shape[2]
    n_replicas=liq_cluster_label.shape[3]
    n_proteins=liq_cluster_size.shape[0]
    
    gyr_radius=np.zeros((n_proteins,n_snapshots,n_replicas))
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            n_clusters=len(np.unique(liq_cluster_label[:,:,i_t,i_replica]))-1
            for i_label in range(n_clusters):
                now_mask=liq_cluster_label[:,:,i_t,i_replica]==i_label+1
                x_ind,y_ind=np.nonzero(now_mask)
                now_rg=0
                for ii in range(len(x_ind)):
                    for jj in range(len(y_ind)):
                        dx = box_size/2 - abs(abs(x_ind[ii]-x_ind[jj]) - box_size/2)
                        dy = box_size/2 - abs(abs(y_ind[ii]-y_ind[jj]) - box_size/2)
                        now_rg+=(dx**2+dy**2)
                now_rg=np.sqrt(now_rg/(2*len(x_ind)**2))
                gyr_radius[i_label,i_t,i_replica]=now_rg
                 
    return gyr_radius



@njit
def percolation_check(liq_cluster_label):
    n_snapshots=liq_cluster_label.shape[2]
    n_replicas=liq_cluster_label.shape[3]
    cluster_hor_max_size=np.zeros((n_snapshots,n_replicas),dtype='int16')
    cluster_ver_max_size=np.zeros((n_snapshots,n_replicas),dtype='int16')
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            n_clusters=len(np.unique(liq_cluster_label[:,:,i_t,i_replica]))-1
            now_max_hor_size=0
            now_max_ver_size=0
            for i_label in range(n_clusters):
                now_mask=liq_cluster_label[:,:,i_t,i_replica]==i_label+1
                now_hor_size=len(np.unique(np.nonzero(now_mask)[0]))
                now_ver_size=len(np.unique(np.nonzero(now_mask)[1]))
                if now_hor_size>now_max_hor_size:
                    now_max_hor_size=now_hor_size
                if now_ver_size>now_max_ver_size:
                    now_max_ver_size=now_ver_size
            cluster_hor_max_size[i_t,i_replica]=now_max_hor_size
            cluster_ver_max_size[i_t,i_replica]=now_max_ver_size
    return cluster_hor_max_size,cluster_ver_max_size


@njit
def HK_cluster_labeling(all_x_pos,all_y_pos,box_size):
    n_monomers=all_x_pos.shape[0]
    n_proteins=all_x_pos.shape[1]
    n_snapshots=all_x_pos.shape[2]
    n_replicas=all_x_pos.shape[3]
    
    aa_present=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype='bool')
    hw_present=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype='bool')
    
    liq_cluster_label=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype='int32')
    prot_belongs_to_liq_cluster=np.zeros((n_proteins,n_snapshots,n_replicas), dtype='int32')
    liq_cluster_size=np.zeros((n_proteins,n_snapshots,n_replicas), dtype='int32')
    
    sol_cluster_label=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype='int32')
    prot_belongs_to_sol_cluster=np.zeros((n_proteins,n_snapshots,n_replicas), dtype='int32')
    sol_cluster_size=np.zeros((n_proteins,n_snapshots,n_replicas), dtype='int32')
    
    liq_cluster_label_nw=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype='int32')
    prot_belongs_to_liq_cluster_nw=np.zeros((n_proteins*n_monomers,n_snapshots,n_replicas), dtype='int32')
    liq_cluster_size_nw=np.zeros((n_proteins*n_monomers,n_snapshots,n_replicas), dtype='int32')
    
    sol_cluster_label_nw=np.zeros((box_size,box_size,n_snapshots,n_replicas), dtype='int32')
    prot_belongs_to_sol_cluster_nw=np.zeros((n_proteins*n_monomers,n_snapshots,n_replicas), dtype='int32')
    sol_cluster_size_nw=np.zeros((n_proteins*n_monomers,n_snapshots,n_replicas), dtype='int32')
    
    
    
    
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            aa_belongs_to=np.zeros((box_size,box_size), dtype='int32')
            for i_prot in range(n_proteins):
                for i_mon in range(n_monomers):
                    aa_present[
                        all_x_pos[i_mon,i_prot,i_t,i_replica]-1,
                        all_y_pos[i_mon,i_prot,i_t,i_replica]-1,i_t,i_replica]=True
                    aa_belongs_to[
                        all_x_pos[i_mon,i_prot,i_t,i_replica]-1,
                        all_y_pos[i_mon,i_prot,i_t,i_replica]-1]=i_prot
                    now_center_x=all_x_pos[i_mon,i_prot,i_t,i_replica]-1
                    now_center_y=all_y_pos[i_mon,i_prot,i_t,i_replica]-1
                    for i_x,i_y in [[-1,0],[1,0],[0,-1],[0,1]]:
                        now_x_ind=((now_center_x+i_x+box_size)%box_size)
                        now_y_ind=((now_center_y+i_y+box_size)%box_size)
                        hw_present[now_x_ind,now_y_ind,i_t,i_replica]=True
            for i_prot in range(n_proteins):
                for i_mon in range(n_monomers):
                    hw_present[all_x_pos[i_mon,i_prot,i_t,i_replica]-1,
                        all_y_pos[i_mon,i_prot,i_t,i_replica]-1,i_t,i_replica]=False

            
            # liquid cluster labelling
            drop_present=np.logical_or(hw_present[:,:,i_t,i_replica],aa_present[:,:,i_t,i_replica])
            
            largest_label=0
            label = np.zeros((box_size,box_size),dtype='int32')
            for i_row in range(box_size):
                for i_col in range(box_size):
                    if drop_present[i_row,i_col]:
                        i_row_u=(i_row-1+box_size)%box_size
                        i_col_l=(i_col-1+box_size)%box_size
                        if (label[i_row_u,i_col] == 0) and (label[i_row,i_col_l] == 0):
                            largest_label+=1
                            label[i_row,i_col]=largest_label
                        elif (label[i_row_u,i_col] != 0) and (label[i_row,i_col_l] == 0):
                            label[i_row,i_col]=label[i_row_u,i_col]
                        elif (label[i_row_u,i_col] == 0) and (label[i_row,i_col_l] != 0):
                            label[i_row,i_col]=label[i_row,i_col_l]
                        else:
                            min_label=min(label[i_row_u,i_col],label[i_row,i_col_l])
                            max_label=max(label[i_row_u,i_col],label[i_row,i_col_l])
                            now_ind_to_change=np.argwhere(label==max_label)
                            for ii in range(now_ind_to_change.shape[0]):
                                label[now_ind_to_change[ii,0],now_ind_to_change[ii,1]]=min_label
                            label[i_row,i_col]=label[i_row,i_col_l]
                            
            now_labels=np.unique(label)
            n_clusters=len(now_labels)-1
            for i_label in range(1,n_clusters+1):
                now_label=now_labels[i_label]
                if i_label != now_label:
                    now_ind_to_change=np.argwhere(label==now_label)
                    for ii in range(now_ind_to_change.shape[0]):
                        label[now_ind_to_change[ii,0],now_ind_to_change[ii,1]]=i_label
                
            
            liq_cluster_label_nw[:,:,i_t,i_replica]=label
            for i_prot in range(n_proteins):
                prot_belongs_to_liq_cluster_nw[i_prot,i_t,i_replica]=label[
                    all_x_pos[0,i_prot,i_t,i_replica]-1,
                    all_y_pos[0,i_prot,i_t,i_replica]-1]
                liq_cluster_size_nw[label[
                    all_x_pos[0,i_prot,i_t,i_replica]-1,
                    all_y_pos[0,i_prot,i_t,i_replica]-1]-1,i_t,i_replica]+=1
                          
            # apply PBCs:
            for i_row in range(box_size):
                if (label[i_row,0] !=0) and (label[i_row,-1] !=0):
                    if (label[i_row,0] != label[i_row,-1]):
                        min_label=min(label[i_row,0],label[i_row,-1])
                        max_label=max(label[i_row,0],label[i_row,-1])
                        now_ind_to_change=np.argwhere(label==max_label)
                        for ii in range(now_ind_to_change.shape[0]):
                            label[now_ind_to_change[ii,0],now_ind_to_change[ii,1]]=min_label
            for i_col in range(box_size):
                if (label[0,i_col] !=0) and (label[-1,i_col] !=0):
                    if (label[0,i_col] != label[-1,i_col]):
                        min_label=min(label[0,i_col],label[-1,i_col])
                        max_label=max(label[0,i_col],label[-1,i_col])
                        now_ind_to_change=np.argwhere(label==max_label)
                        for ii in range(now_ind_to_change.shape[0]):
                            label[now_ind_to_change[ii,0],now_ind_to_change[ii,1]]=min_label
                            
            now_labels=np.unique(label)
            n_clusters=len(now_labels)-1
            for i_label in range(1,n_clusters+1):
                now_label=now_labels[i_label]
                if i_label != now_label:
                    now_ind_to_change=np.argwhere(label==now_label)
                    for ii in range(now_ind_to_change.shape[0]):
                        label[now_ind_to_change[ii,0],now_ind_to_change[ii,1]]=i_label
                
            
            liq_cluster_label[:,:,i_t,i_replica]=label
            for i_prot in range(n_proteins):
                prot_belongs_to_liq_cluster[i_prot,i_t,i_replica]=label[
                    all_x_pos[0,i_prot,i_t,i_replica]-1,
                    all_y_pos[0,i_prot,i_t,i_replica]-1]
                liq_cluster_size[label[
                    all_x_pos[0,i_prot,i_t,i_replica]-1,
                    all_y_pos[0,i_prot,i_t,i_replica]-1]-1,i_t,i_replica]+=1
                
            # solid cluster labelling
            largest_label=0
            label = np.zeros((box_size,box_size),dtype='int32')
            for i_row in range(box_size):
                for i_col in range(box_size):
                    if aa_present[i_row,i_col,i_t,i_replica]:
                        i_row_u=(i_row-1+box_size)%box_size
                        i_col_l=(i_col-1+box_size)%box_size
                        if (label[i_row_u,i_col] == 0) and (label[i_row,i_col_l] == 0):
                            largest_label+=1
                            label[i_row,i_col]=largest_label
                        elif (label[i_row_u,i_col] != 0) and (label[i_row,i_col_l] == 0):
                            label[i_row,i_col]=label[i_row_u,i_col]
                        elif (label[i_row_u,i_col] == 0) and (label[i_row,i_col_l] != 0):
                            label[i_row,i_col]=label[i_row,i_col_l]
                        else:
                            min_label=min(label[i_row_u,i_col],label[i_row,i_col_l])
                            max_label=max(label[i_row_u,i_col],label[i_row,i_col_l])
                            now_ind_to_change=np.argwhere(label==max_label)
                            for ii in range(now_ind_to_change.shape[0]):
                                label[now_ind_to_change[ii,0],now_ind_to_change[ii,1]]=min_label
                            label[i_row,i_col]=label[i_row,i_col_l]
                            
            now_labels=np.unique(label)
            n_clusters=len(now_labels)-1
            for i_label in range(1,n_clusters+1):
                now_label=now_labels[i_label]
                if i_label != now_label:
                    now_ind_to_change=np.argwhere(label==now_label)
                    for ii in range(now_ind_to_change.shape[0]):
                        label[now_ind_to_change[ii,0],now_ind_to_change[ii,1]]=i_label
            
            sol_cluster_label_nw[:,:,i_t,i_replica]=label
            
            for i_prot in range(n_proteins):
                prot_belongs_to_sol_cluster_nw[i_prot,i_t,i_replica]=label[
                    all_x_pos[0,i_prot,i_t,i_replica]-1,
                    all_y_pos[0,i_prot,i_t,i_replica]-1]
                sol_cluster_size_nw[label[
                    all_x_pos[0,i_prot,i_t,i_replica]-1,
                    all_y_pos[0,i_prot,i_t,i_replica]-1]-1,i_t,i_replica]+=1
                
            #apply PBCs:
            for i_row in range(box_size):
                if (label[i_row,0] !=0) and (label[i_row,-1] !=0):
                    if (label[i_row,0] != label[i_row,-1]):
                        min_label=min(label[i_row,0],label[i_row,-1])
                        max_label=max(label[i_row,0],label[i_row,-1])
                        now_ind_to_change=np.argwhere(label==max_label)
                        for ii in range(now_ind_to_change.shape[0]):
                            label[now_ind_to_change[ii,0],now_ind_to_change[ii,1]]=min_label
        
            for i_col in range(box_size):
                if (label[0,i_col] !=0) and (label[-1,i_col] !=0):
                    if (label[0,i_col] != label[-1,i_col]):
                        min_label=min(label[0,i_col],label[-1,i_col])
                        max_label=max(label[0,i_col],label[-1,i_col])
                        now_ind_to_change=np.argwhere(label==max_label)
                        for ii in range(now_ind_to_change.shape[0]):
                            label[now_ind_to_change[ii,0],now_ind_to_change[ii,1]]=min_label
                            
            now_labels=np.unique(label)
            n_clusters=len(now_labels)-1
            for i_label in range(1,n_clusters+1):
                now_label=now_labels[i_label]
                if i_label != now_label:
                    now_ind_to_change=np.argwhere(label==now_label)
                    for ii in range(now_ind_to_change.shape[0]):
                        label[now_ind_to_change[ii,0],now_ind_to_change[ii,1]]=i_label
            
            sol_cluster_label[:,:,i_t,i_replica]=label
            
            for i_prot in range(n_proteins):
                prot_belongs_to_sol_cluster[i_prot,i_t,i_replica]=label[
                    all_x_pos[0,i_prot,i_t,i_replica]-1,
                    all_y_pos[0,i_prot,i_t,i_replica]-1]
                sol_cluster_size[label[
                    all_x_pos[0,i_prot,i_t,i_replica]-1,
                    all_y_pos[0,i_prot,i_t,i_replica]-1]-1,i_t,i_replica]+=1
            
    return liq_cluster_label,liq_cluster_size,prot_belongs_to_liq_cluster,\
        sol_cluster_label,sol_cluster_size,prot_belongs_to_sol_cluster,\
        liq_cluster_label_nw,liq_cluster_size_nw,prot_belongs_to_liq_cluster_nw,\
        sol_cluster_label_nw,sol_cluster_size_nw,prot_belongs_to_sol_cluster_nw
        
@njit
def from_cluster_sizes_to_cluster_volumes(
    liq_cluster_label,liq_cluster_size,
    sol_cluster_label,sol_cluster_size,
    V_cell):
    n_snapshots=liq_cluster_label.shape[2]
    n_replicas=liq_cluster_label.shape[3]
    
    tot_liq_cluster_sizes=np.zeros(liq_cluster_size.shape,dtype='int32')
    tot_liq_cluster_volumes=np.zeros(liq_cluster_size.shape)
    tot_liq_cluster_volume_fractions=np.zeros(liq_cluster_size.shape)
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            now_n_clusters=len(np.nonzero(liq_cluster_size[:,i_t,i_replica])[0])
            for ii in range(now_n_clusters):
                now_label_mask=liq_cluster_label[:,:,i_t,i_replica]==(ii+1)
                tot_liq_cluster_sizes[ii,i_t,i_replica]=now_label_mask.sum()
                for ix in range(now_label_mask.shape[0]):
                    for iy in range(now_label_mask.shape[1]):
                        if now_label_mask[ix,iy]:
                            tot_liq_cluster_volumes[ii,i_t,i_replica]+=V_cell[ix,iy,i_t,i_replica]
                tot_liq_cluster_volume_fractions[ii,i_t,i_replica]=tot_liq_cluster_volumes[ii,i_t,i_replica]/(V_cell[:,:,i_t,i_replica].sum())
    
    
    
    tot_sol_cluster_sizes=np.zeros(sol_cluster_size.shape,dtype='int32')
    tot_sol_cluster_volumes=np.zeros(liq_cluster_size.shape)
    tot_sol_cluster_volume_fractions=np.zeros(liq_cluster_size.shape)
    for i_replica in range(n_replicas):
        for i_t in range(n_snapshots):
            now_n_clusters=len(np.nonzero(sol_cluster_size[:,i_t,i_replica])[0])
            for ii in range(now_n_clusters):
                now_label_mask=sol_cluster_label[:,:,i_t,i_replica]==(ii+1)
                tot_sol_cluster_sizes[ii,i_t,i_replica]=now_label_mask.sum()
                for ix in range(now_label_mask.shape[0]):
                    for iy in range(now_label_mask.shape[1]):
                        if now_label_mask[ix,iy]:
                            tot_sol_cluster_volumes[ii,i_t,i_replica]+=V_cell[ix,iy,i_t,i_replica]
                tot_sol_cluster_volume_fractions[ii,i_t,i_replica]=tot_sol_cluster_volumes[ii,i_t,i_replica]/(V_cell[:,:,i_t,i_replica].sum())
    return tot_liq_cluster_sizes,tot_liq_cluster_volumes,tot_liq_cluster_volume_fractions,\
            tot_sol_cluster_sizes,tot_sol_cluster_volumes,tot_sol_cluster_volume_fractions



