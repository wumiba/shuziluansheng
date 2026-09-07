# -*- coding: utf-8 -*-
import sys, os, struct
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import _fonts

def load_tri(filename):
    with open(filename,'rb') as f:
        hdr = f.read(80)
        n = struct.unpack('<I', f.read(4))[0]
        buf = np.frombuffer(f.read(n*50), dtype=np.uint8).reshape(n, 50)
    a = buf[:, :48].view('<f4').reshape(n, 12)   # 3 normal + 9 vertex (跳过2B attr)
    tri = a[:, 3:12].reshape(n,3,3).copy()
    return tri

chest = load_tri('d.stl')   # 胸腔
heart = load_tri('b.stl')   # 心脏
step = max(1, len(heart)//8000)
heart = heart[::step]
print('chest faces', len(chest), 'heart faces (decimated)', len(heart))

def normals(tri):
    v0,v1,v2 = tri[:,0],tri[:,1],tri[:,2]
    n = np.cross(v1-v0, v2-v0)
    nl = np.linalg.norm(n, axis=1, keepdims=True)+1e-12
    return n/nl

L = np.array([0.4,0.9,0.5]); L/=np.linalg.norm(L)

def shade(base_rgb, tri):
    n = normals(tri)
    k = 0.6+0.4*np.clip(n@L, 0, 1)      # 平坦光照
    return (np.array(base_rgb)[None,:]*k[:,None])

def render(ax, tri, rgb, alpha, crop=None, edge=False, lw=0.1):
    if crop is not None:
        c = tri.mean(axis=1)
        m = (np.abs(c[:,0])<=crop[0]) & (np.abs(c[:,1])<=crop[1]) & (np.abs(c[:,2])<=crop[2])
        tri = tri[m]
    cols = shade(rgb, tri)
    coll = Poly3DCollection(tri, facecolors=np.concatenate([cols, np.full((len(tri),1), alpha)], axis=1),
                            edgecolors=(0.2,0.2,0.2,0.4) if edge else 'none', linewidths=lw, zsort='average')
    ax.add_collection3d(coll)

def set_eq(ax, tri_list):
    los = np.array([t.reshape(-1,3).min(0) for t in tri_list])
    his = np.array([t.reshape(-1,3).max(0) for t in tri_list])
    lo, hi = los.min(0), his.max(0)
    c = (lo+hi)/2; r = float(np.max(hi-lo))/2
    ax.set_xlim(c[0]-r, c[0]+r); ax.set_ylim(c[1]-r, c[1]+r); ax.set_zlim(c[2]-r, c[2]+r)

fig = plt.figure(figsize=(12.4, 5.6))

# ---- (a) 全局视图 ----
ax = fig.add_subplot(121, projection='3d')
ax.set_axis_off()
render(ax, chest, (0.86,0.88,0.93), 0.30)
render(ax, heart, (0.84,0.16,0.18), 0.99)
set_eq(ax, [chest, heart])
ax.view_init(elev=14, azim=-65)
ax.text2D(0.02, 0.94, '(a) 胸腔与心脏三维几何模型（同一坐标系，胸腔半透明显示）',
          transform=ax.transAxes, fontsize=10.5, color='#111')

# ---- (b) 中心剖视（胸腔壁切片内可见心脏）----
ax2 = fig.add_subplot(122, projection='3d')
ax2.set_axis_off()
crop = (0.055, 0.09, 0.10)   # 只保留胸壁中心区域切片
render(ax2, chest, (0.86,0.88,0.93), 0.34, crop=crop)
render(ax2, heart, (0.84,0.16,0.18), 1.0)
ax2.view_init(elev=14, azim=-65)
ax2.text2D(0.02, 0.94, '(b) 中心剖视：胸腔壁切片包裹心脏（几何孪生对象）',
          transform=ax2.transAxes, fontsize=10.5, color='#111')

fig.savefig('数字孪生素材_开题报告/图/fig2_胸腔心脏三维几何模型.png', dpi=300)
print('saved fig2')
