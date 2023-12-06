from analysis_tools.animation import Animator

path = 'simulation_results/jupiterlike_planet/'
animator = Animator(path, xlabel='x', ylabel='y', xlim=0.01, ylim=0.01)
animator.make_animation(save_path='simulation_results/animation_jup.mp4')
