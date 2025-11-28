import os
import json
import requests
import base64
from typing import Dict, Any, Optional, Union, List
from PIL import Image
import io

def encode_image_to_base64(image_path: str) -> str:
    """
    Encode an image file to base64 string.
    
    Args:
        image_path: Path to the image file
        
    Returns:
        str: Base64 encoded image string
    """
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode('utf-8')

def get_image_mime_type(image_path: str) -> str:
    """
    Get the MIME type of an image file.
    
    Args:
        image_path: Path to the image file
        
    Returns:
        str: MIME type (e.g., 'image/jpeg', 'image/png')
    """
    ext = os.path.splitext(image_path)[1].lower()
    mime_types = {
        '.jpg': 'image/jpeg',
        '.jpeg': 'image/jpeg',
        '.png': 'image/png',
        '.gif': 'image/gif',
        '.bmp': 'image/bmp',
        '.webp': 'image/webp'
    }
    return mime_types.get(ext, 'image/jpeg')

def call_llm_with_image(
    prompt: str,
    image_input: Union[str, List[str]],
    provider: str = "openai",
    model: str = None,
    api_key: Optional[str] = None,
    temperature: float = 0.7,
    max_tokens: int = 4096,
    system_prompt: Optional[str] = None,
    additional_params: Optional[Dict[str, Any]] = None
) -> str:
    """
    Call a vision-enabled LLM from various providers with image input.
    
    Args:
        prompt: The user prompt to send to the LLM
        image_input: Single image path/URL or list of image paths/URLs
        provider: One of "openai", "anthropic", or "openrouter", or custom URL
        model: Specific vision model from the provider to use
        api_key: API key for the provider (if None, gets from environment)
        temperature: Sampling temperature (0-1)
        max_tokens: Maximum tokens to generate
        system_prompt: Optional system prompt for providers that support it
        additional_params: Additional parameters to pass to the provider's API
    
    Returns:
        str: The generated text response
    """
    provider = provider.lower()
    additional_params = additional_params or {}
    
    # Default vision models for each provider if not specified
    default_models = {
        "openai": "gpt-4o",
        "anthropic": "claude-3-5-sonnet-20241022",
        "openrouter": "openai/gpt-4o",
    }
    
    # Use default model if not specified
    if not model:
        model = default_models.get(provider)
        if not model:
            raise ValueError(f"No model specified and no default available for provider: {provider}")
    
    # Get API key from environment if not provided
    if not api_key:
        env_var_names = {
            "openai": "OPENAI_API_KEY",
            "anthropic": "ANTHROPIC_API_KEY",
            "openrouter": "OPENROUTER_API_KEY",
        }
        env_var = env_var_names.get(provider)
        if env_var:
            api_key = os.environ.get(env_var)
            if not api_key:
                raise ValueError(f"API key not provided and {env_var} not found in environment")
    
    # Ensure image_input is a list
    if isinstance(image_input, str):
        image_input = [image_input]
    
    # Prepare messages format
    messages = []
    if system_prompt:
        messages.append({"role": "system", "content": system_prompt})
    
    # OpenAI API call
    if provider == "openai":
        try:
            import openai
        except ImportError:
            raise ImportError("Please install openai package: pip install openai")
        
        client = openai.OpenAI(api_key=api_key)
        
        # Prepare content with images
        content = [{"type": "text", "text": prompt}]
        
        for image_path in image_input:
            if image_path.startswith("http"):
                # URL image
                content.append({
                    "type": "image_url",
                    "image_url": {"url": image_path}
                })
            else:
                # Local image file
                if not os.path.exists(image_path):
                    raise FileNotFoundError(f"Image file not found: {image_path}")
                
                base64_image = encode_image_to_base64(image_path)
                mime_type = get_image_mime_type(image_path)
                
                content.append({
                    "type": "image_url",
                    "image_url": {
                        "url": f"data:{mime_type};base64,{base64_image}"
                    }
                })
        
        messages.append({"role": "user", "content": content})
        
        response = client.chat.completions.create(
            model=model,
            messages=messages,
            temperature=temperature,
            max_tokens=max_tokens,
            **additional_params
        )
        
        return response.choices[0].message.content
    
    # Custom OpenAI-compatible API call (base_url as provider)
    elif provider.startswith("http"):
        try:
            import openai
        except ImportError:
            raise ImportError("Please install openai package: pip install openai")
        
        custom_api_key = api_key or os.environ.get("CUSTOMIZED_API_KEY")
        if not custom_api_key:
            raise ValueError("API key not provided and CUSTOMIZED_API_KEY not found in environment")
        
        client = openai.OpenAI(api_key=custom_api_key, base_url=provider)
        
        # Prepare content with images
        content = [{"type": "text", "text": prompt}]
        
        for image_path in image_input:
            if image_path.startswith("http"):
                # URL image
                content.append({
                    "type": "image_url",
                    "image_url": {"url": image_path}
                })
            else:
                # Local image file
                if not os.path.exists(image_path):
                    raise FileNotFoundError(f"Image file not found: {image_path}")
                
                base64_image = encode_image_to_base64(image_path)
                mime_type = get_image_mime_type(image_path)
                
                content.append({
                    "type": "image_url",
                    "image_url": {
                        "url": f"data:{mime_type};base64,{base64_image}"
                    }
                })
        
        # Handle message history properly
        api_messages = messages.copy()
        api_messages.append({"role": "user", "content": content})
        
        # If additional_params contains message history, merge it properly
        if 'messages' in additional_params:
            # Use the full conversation history from additional_params instead
            history_messages = additional_params.pop('messages')
            api_messages = history_messages
            
            # Add the current message with image to the history
            api_messages.append({"role": "user", "content": content})
            
            # Only add system prompt if it's not already in the history
            if system_prompt and not any(msg.get('role') == 'system' for msg in api_messages):
                api_messages.insert(0, {"role": "system", "content": system_prompt})
        
        # Call the API with the proper message history
        response = client.chat.completions.create(
            model=model,
            messages=api_messages,
            temperature=temperature,
            max_tokens=max_tokens,
            **additional_params
        )
        return response.choices[0].message.content
    
    # Anthropic API call
    elif provider == "anthropic":
        try:
            import anthropic
        except ImportError:
            raise ImportError("Please install anthropic package: pip install anthropic")
        
        client = anthropic.Anthropic(api_key=api_key)
        
        # Prepare content with images for Anthropic
        content = [{"type": "text", "text": prompt}]
        
        for image_path in image_input:
            if image_path.startswith("http"):
                # For URLs, we need to download and encode
                response = requests.get(image_path)
                response.raise_for_status()
                base64_image = base64.b64encode(response.content).decode('utf-8')
                
                # Try to determine mime type from headers or URL
                content_type = response.headers.get('content-type', 'image/jpeg')
                if not content_type.startswith('image/'):
                    content_type = 'image/jpeg'
            else:
                # Local image file
                if not os.path.exists(image_path):
                    raise FileNotFoundError(f"Image file not found: {image_path}")
                
                base64_image = encode_image_to_base64(image_path)
                content_type = get_image_mime_type(image_path)
            
            content.append({
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": content_type,
                    "data": base64_image
                }
            })
        
        # Create the message with system as a string
        message_params = {
            "model": model,
            "max_tokens": max_tokens,
            "temperature": temperature,
            "messages": [
                {
                    "role": "user", 
                    "content": content
                }
            ]
        }
        
        # Add system prompt if provided
        if system_prompt:
            message_params["system"] = system_prompt
            
        # Add any additional parameters
        message_params.update(additional_params)
        
        # Call the API
        response = client.messages.create(**message_params)
        
        # Extract the text content from the response
        if hasattr(response, 'content') and len(response.content) > 0:
            content_block = response.content[0]
            if hasattr(content_block, 'text'):
                return content_block.text
            elif isinstance(content_block, dict) and 'text' in content_block:
                return content_block['text']
            else:
                return str(response.content)
        else:
            return "No content returned from Anthropic API"
    
    # OpenRouter API call
    elif provider == "openrouter":
        url = "https://openrouter.ai/api/v1/chat/completions"
        
        headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        }
        
        # Prepare content with images
        content = [{"type": "text", "text": prompt}]
        
        for image_path in image_input:
            if image_path.startswith("http"):
                # URL image
                content.append({
                    "type": "image_url",
                    "image_url": {"url": image_path}
                })
            else:
                # Local image file
                if not os.path.exists(image_path):
                    raise FileNotFoundError(f"Image file not found: {image_path}")
                
                base64_image = encode_image_to_base64(image_path)
                mime_type = get_image_mime_type(image_path)
                
                content.append({
                    "type": "image_url",
                    "image_url": {
                        "url": f"data:{mime_type};base64,{base64_image}"
                    }
                })
        
        messages.append({"role": "user", "content": content})
        
        data = {
            "model": model,
            "messages": messages,
            "temperature": temperature,
            "max_tokens": max_tokens,
            **additional_params
        }
        
        response = requests.post(url, headers=headers, data=json.dumps(data))
        response.raise_for_status()
        
        return response.json()["choices"][0]["message"]["content"]
    
    else:
        raise ValueError(f"Unsupported provider: {provider}")

def call_llm_analyze_image(
    image_input: Union[str, List[str]],
    analysis_prompt: str = "Please analyze this image and describe what you see.",
    provider: str = "openai",
    model: str = None,
    api_key: Optional[str] = None,
    temperature: float = 0.7,
    max_tokens: int = 4096,
    system_prompt: Optional[str] = None,
    additional_params: Optional[Dict[str, Any]] = None
) -> str:
    """
    Convenient function to analyze images with a default prompt.
    
    Args:
        image_input: Single image path/URL or list of image paths/URLs
        analysis_prompt: The prompt for analyzing the image
        provider: One of "openai", "anthropic", or "openrouter", or custom URL
        model: Specific vision model from the provider to use
        api_key: API key for the provider (if None, gets from environment)
        temperature: Sampling temperature (0-1)
        max_tokens: Maximum tokens to generate
        system_prompt: Optional system prompt for providers that support it
        additional_params: Additional parameters to pass to the provider's API
    
    Returns:
        str: The generated analysis response
    """
    return call_llm_with_image(
        prompt=analysis_prompt,
        image_input=image_input,
        provider=provider,
        model=model,
        api_key=api_key,
        temperature=temperature,
        max_tokens=max_tokens,
        system_prompt=system_prompt,
        additional_params=additional_params
    )

def call_llm_extract_text_from_image(
    image_input: Union[str, List[str]],
    provider: str = "openai",
    model: str = None,
    api_key: Optional[str] = None,
    temperature: float = 0.1,
    max_tokens: int = 4096,
    additional_params: Optional[Dict[str, Any]] = None
) -> str:
    """
    Convenient function to extract text from images (OCR).
    
    Args:
        image_input: Single image path/URL or list of image paths/URLs
        provider: One of "openai", "anthropic", or "openrouter", or custom URL
        model: Specific vision model from the provider to use
        api_key: API key for the provider (if None, gets from environment)
        temperature: Sampling temperature (0-1), set low for OCR
        max_tokens: Maximum tokens to generate
        additional_params: Additional parameters to pass to the provider's API
    
    Returns:
        str: The extracted text from the image(s)
    """
    ocr_prompt = """Extract all text from this image. 
    Please provide the text exactly as it appears, preserving formatting and structure where possible.
    If there are multiple columns or sections, please maintain their organization.
    Only return the text content, no additional commentary."""
    
    return call_llm_with_image(
        prompt=ocr_prompt,
        image_input=image_input,
        provider=provider,
        model=model,
        api_key=api_key,
        temperature=temperature,
        max_tokens=max_tokens,
        system_prompt="You are an expert OCR assistant. Extract text accurately and preserve formatting.",
        additional_params=additional_params
    ) 