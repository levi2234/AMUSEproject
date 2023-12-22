from analysis_tools.animation import Animator

path = 'simulation_results/gas_results/4.0e+42_erg/'
save_path = 'simulation_results/animation_4.0e+42_erg.mp4'
animator = Animator(path, xlabel='x', ylabel='y', center='planet', xlim=0.02, ylim=0.02)
animator.make_animation(save_path=save_path)
