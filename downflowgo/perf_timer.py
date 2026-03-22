def runtime(start_time, end_time):
    # Calculate the duration
    execution_time = end_time - start_time

    # Format the execution time
    if execution_time >= 60:
        minutes = int(execution_time // 60)
        seconds = int(execution_time % 60)
        print(f"The code took {minutes} minutes and {seconds} seconds to execute.")
    else:
        print(f"The code took {int(execution_time)} seconds to execute.")