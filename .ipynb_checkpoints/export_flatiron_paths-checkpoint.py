from one.api import ONE
one = ONE()
import numpy as np
from pathlib import Path
from brainbox.io.one import SpikeSortingLoader
import pickle
import requests
import json
import csv
from ibllib.atlas import AllenAtlas
atlas = AllenAtlas(25)
import pandas as pd

uuid_2_index_fp = './data/uuid_index.pkl'
bregma = [5739., 5400.,  332.]

token = one.alyx._par.TOKEN[one.alyx._par.ALYX_LOGIN]

headers = {
    'Authorization': f'Token {list(token.values())[0]}',
    'Accept': 'application/json',
    'Content-Type': 'application/json'}

def export_paths(eid):
    save_path = Path('./data/flatiron_paths')
    probes = ['probe00','probe01']
    alyx_url = one.alyx.base_url
    
    probe_datas = []
    
    insertions = one.alyx.rest('insertions', 'list', session=eid)

    pids = [i['id'] for i in insertions]
    for pid in pids:

        _, probe = one.pid2eid(pid)
        if probe==probes[0]:
            pi = 0
        else:
            pi = 1

        # List of datasets, first column is dataset type, second column is collection
        dsets = [('spikes.times', f'alf/{probe}'),
                 ('spikes.clusters', f'alf/{probe}'),
                 ('clusters.mlapdv', f'alf/{probe}'),  # This is only available for resolved sessions
                 ('wheel.position', 'alf'),
                 ('wheel.timestamps', 'alf'),
                 ('_iblrig_stimPositionScreen.raw', 'raw_behavior_data'), # I think this is the correct file but i'm not sure it
                 ('trials.goCue_times', 'alf'), # time of tone to indicate trial start
                 ('trials.feedback_times', 'alf'), # time of valve click or white noise to indicate reward/ punishment
                 ('trials.feedbackType', 'alf'), # whether correct or incorrect
                 ('trials.contrastLeft', 'alf'),
                 ('trials.contrastRight', 'alf')
        ]
        
        # make sure we at a minimum have spikes
        dataset = one.alyx.rest('datasets', 'list', session=eid,
            dataset_type=['spikes.times'], collection=f'alf/{probe}')
        if len(dataset)>0:

            urls = []
            for dset in dsets:
                dataset = one.alyx.rest('datasets', 'list', session=eid,
                                        dataset_type=dset[0], collection=dset[1])
                if len(dataset) > 0:
                    urls.append(next(fr['data_url'] for fr in dataset[0]['file_records']
                                     if 'flatiron' in fr['data_repository']))

            outfile = save_path.joinpath(f'file_urls_{eid}_{probe}.txt')

            with open(outfile, 'w') as f:
                f.write(eid)
                f.write('\n')
                f.write(str(pi))
                f.write('\n')
                f.write(pid)
                f.write('\n')
                for url in urls:
                    f.write(url)
                    f.write('\n')
                f.close()
                        
            # get the mlapdv coordinates and save those
            ssl = SpikeSortingLoader(pid=pid, one=one)
            _, clusters, channels = ssl.load_spike_sorting()
            print(f"{pid} {ssl.pname} {ssl.collection} {ssl.histology}")
            
            # try converting channel positions to mlapdv
            x = channels['x']
            y = channels['y']
            z = channels['z']
            xyz_coords = []
            for xv,yv,zv in zip(x,y,z):
                xyz_coords.append([xv,yv,zv])
            mlapdv = atlas.xyz2ccf(xyz_coords)
            clusters_coords = mlapdv[clusters.channels]
            # save coordinates to a csv file
            pd.DataFrame(clusters_coords,columns=['ml','ap','dv']).to_csv('data/'+pid+'.csv')
            
            # also load the probe trajectory data and save that to a file as well
            print(eid + " " + probe)
            query = f'/trajectories?session={eid}&probe={probe}&provenance=Ephys+aligned+histology+track'
            r = requests.get(alyx_url + query, headers=headers)
            traj = json.loads(r.text)
            res = traj['results']
            if len(res)>0:
                traj_data = traj['results'][0]
                xyz = [traj_data['x']/1000000, traj_data['y']/1000000, traj_data['z']/1000000]
                mlapdv = atlas.xyz2ccf(xyz)
                probe_data = [eid, probe, mlapdv[0], mlapdv[1], mlapdv[2],  traj_data['depth'], traj_data['theta'], traj_data['phi'], traj_data['roll']]
                probe_datas.append(probe_data)
            else:
                print("MISSING TRAJECTORY")
    return probe_datas
        
def export_paths_sessions(sess_id):
    adatas = []
    for eid in sess_id:
        probe_datas = export_paths(eid)
        for probe_data in probe_datas:
            adatas.append(probe_data)
            
    fields = ['eid','probe','ml','ap','dv','depth','theta','phi','roll']
    with open('./data/flatiron_paths/probe_trajectories.csv', 'w') as f:

        # using csv.writer method from CSV package
        write = csv.writer(f)

        write.writerow(fields)
        write.writerows(adatas)
        

# testing

sess_id = ['73918ae1-e4fd-4c18-b132-00cb555b1ad2']
# Run
export_paths_sessions(sess_id)