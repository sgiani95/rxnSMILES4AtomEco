def get_temperature(lower_limit=-273.15, upper_limit=1000):
    try:
        # Request temperature from user
        temp_C = float(input("Please enter the temperature in Celsius (°C): "))
        
        # Validate that the temperature is within the specified limits
        if temp_C < lower_limit:
            print(f"Temperature cannot be below {lower_limit}°C.")
            return None
        if temp_C > upper_limit:
            print(f"Temperature cannot exceed {upper_limit}°C.")
            return None
        
        return temp_C
    except ValueError:
        # Handle the case where the input is not a valid number
        print("Invalid input. Please enter a valid number.")
        return None

# Example usage:
temperature = get_temperature()
if temperature is not None:
    print(f"The temperature you entered is {temperature}°C.")
