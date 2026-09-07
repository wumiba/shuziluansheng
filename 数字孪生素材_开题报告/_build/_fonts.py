from matplotlib import font_manager, rcParams
for fp in [r'C:\Windows\Fonts\simhei.ttf', r'C:\Windows\Fonts\msyh.ttc']:
    try:
        font_manager.fontManager.addfont(fp)
    except Exception as e:
        print('skip', fp, e)
names = [f.name for f in font_manager.fontManager.ttflist if f.name in ('SimHei','Microsoft YaHei')]
print('available CJK fonts:', sorted(set(names)))
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Microsoft YaHei','SimHei']
rcParams['axes.unicode_minus'] = False
