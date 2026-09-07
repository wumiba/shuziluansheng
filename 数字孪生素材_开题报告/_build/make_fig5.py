# -*- coding: utf-8 -*-
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PIL import Image
import numpy as np
import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import _fonts

im = Image.open('chest_heart.gif')
n = im.n_frames
idx = [0, n//3, (2*n)//3] if n >= 3 else list(range(n))
frames = []
for i in idx:
    im.seek(i); frames.append(im.convert('RGB'))

fig, axs = plt.subplots(1, len(frames), figsize=(6.2*len(frames), 5.0))
if len(frames)==1: axs=[axs]
for ax, fr, i in zip(axs, frames, idx):
    ax.imshow(fr)
    ax.axis('off')
    ax.set_title(f'(仿真动画帧 · 第 {i+1} 帧)', fontsize=11)
fig.suptitle('胸腔呼吸＋心脏搏动联合形变动画（数字孪生仿真输出，选自随附 GIF）', fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig('数字孪生素材_开题报告/图/fig5_胸腔心脏仿真动画帧.png', dpi=200)
# 单帧示例
mid = frames[len(frames)//2]
mid.save('数字孪生素材_开题报告/图/fig5_胸腔心脏仿真动画_单帧示例.png')
print('fig5 done, frames at idx', idx)
