# Facility Booking System

## Problem Statement

There are N computer centers, each with multiple facility rooms offering various types of facilities. Users can request bookings for these rooms across 24 time slots in a day. Each booking specifies the computer center, facility room, starting slot, and number of consecutive slots required. Requests are processed concurrently using GPU cores to maximize throughput and performance.

## System Requirements

- Users interact with an application to book facility rooms:
  - Specify computer center, facility room, starting time slot, and number of slots.
  - Slots are booked consecutively from the starting slot.
  
- **Concurrency Handling**:
  - Ensure consistency and avoid race conditions when multiple requests target the same facility.
  - Prioritize requests with smaller IDs if multiple requests target the same facility and the room is at capacity.
  
- **GPU Utilization**:
  - Implement parallel processing using GPU kernels for efficient handling of multiple requests.
  - Avoid print statements inside GPU kernels to maintain performance.

## Output Requirements

After processing all requests:
- Print the total number of successful and failed requests.
- Print the total successful and failed requests for each computer center.

## Additional Points

- **Race Condition**: Ensure mutual exclusion to prevent race conditions when multiple requests access the same facility concurrently.
- **Priority Handling**: Requests with smaller IDs should be prioritized if the facility is overbooked.
- **GPU Implementation**: Use GPU kernels effectively to achieve parallelism and improve performance.
- **Sequential Implementation**: Minimize sequential execution to maximize marks; focus on leveraging GPU parallelism.

## Submission Guidelines

- Submit a `.cu` file named `RollNumber.cu` on Moodle.
  
## Learning Suggestions

- Develop a CPU version of the code to compare performance with the GPU implementation.
- Measure and analyze performance differences between CPU and GPU implementations using large inputs.
- Explore synchronization mechanisms to manage concurrent access to shared resources effectively.

