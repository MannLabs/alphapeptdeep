if __name__ == "__main__":
    try:
        import peptdeep.gui
        import multiprocessing
        multiprocessing.freeze_support()
        peptdeep.gui.run()
    except KeyboardInterrupt:
        pass
    except Exception:
        import traceback
        import sys
        exc_info = sys.exc_info()
        # Display the *original* exception
        traceback.print_exception(*exc_info)
        input("Something went wrong, press any key to continue...")
